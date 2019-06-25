#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <dirent.h>
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
#include "mclib_3d.h"
#include <omp.h>

#define R_DIM 1260
#define THETA_DIM 280
#define PHI_DIM 280


void read_hydro(char hydro_prefix[200], int frame, double r_inj, double **x, double **y,  double **z, double **szx, double **szy, double **r, double **theta, double **phi,\
  double **velx, double **vely, double **velz, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number,  int ph_inj_switch, double min_r, double max_r, double fps, FILE *fPtr)
{
    FILE *hydroPtr=NULL;
    char hydrofile[200]="", file_num[200]="", full_file[200]="",file_end[200]=""  ;
    char buf[10]="";
    int i=0, j=0, k=0, elem=0, elem_factor=0;
    int phi_min_index=0, phi_max_index=0, r_min_index=0, r_max_index=0, theta_min_index=0, theta_max_index=0; //all_index_buffer contains phi_min, phi_max, theta_min, theta_max, r_min, r_max indexes to get from grid files
    int r_index=0, theta_index=0, phi_index=0, hydro_index=0, all_index_buffer=0, adjusted_remapping_index=0, dr_index=0;
    int *remapping_indexes=NULL;
    float buffer=0;
    float *dens_unprc=NULL;
    float *vel_r_unprc=NULL;
    float *vel_theta_unprc=NULL;
    float *vel_phi_unprc=NULL;
    float *pres_unprc=NULL;
    double ph_rmin=0, ph_rmax=0;
    double r_in=1e10, r_ref=2e13;
    double *r_edge=NULL;
    double *dr=NULL;
    double *r_unprc=malloc(sizeof(double)*R_DIM);
    double *theta_unprc=malloc(sizeof(double)*THETA_DIM);
    double *phi_unprc=malloc(sizeof(double)*PHI_DIM);
    
    if (ph_inj_switch==0)
    {
        ph_rmin=min_r;
        ph_rmax=max_r;
    }
    
    snprintf(file_end,sizeof(file_end),"%s","small.data" );
    
    //density
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 1,"-" );
    modifyFlashName(file_num, hydrofile, frame,1);
    
    fprintf(fPtr,">> Opening file %s\n", file_num);
    fflush(fPtr);
    
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    /*
    fprintf(fPtr,"Reading Density: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(&phi_min_index, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&phi_max_index, sizeof(int)*1, 1,hydroPtr);
    fread(&theta_min_index, sizeof(int)*1, 1,hydroPtr);
    fread(&theta_max_index, sizeof(int)*1, 1,hydroPtr);
    fread(&r_min_index, sizeof(int)*1, 1,hydroPtr);
    fread(&r_max_index, sizeof(int)*1, 1,hydroPtr);
    fclose(hydroPtr);
    
    //fortran indexing starts @ 1, but C starts @ 0
    r_min_index--;
    r_max_index--;
    theta_min_index--;
    theta_max_index--;
    phi_min_index--;
    phi_max_index--;
    
    //number of elements defined by this now
    elem=(r_max_index+1-r_min_index)*(theta_max_index+1-theta_min_index)*(phi_max_index+1-phi_min_index); //add 1 b/c max_index is 1 less than max number of elements in file
    /*
    fprintf(fPtr,"Elem %d\n", elem);
    fprintf(fPtr,"Limits %d, %d, %d, %d, %d, %d\n", phi_min_index, phi_max_index, theta_min_index, theta_max_index, r_min_index, r_max_index); 
    fflush(fPtr);
    */
    //now with number of elements allocate data, remember last element is some garbage that only fortran uses
    dens_unprc=malloc(elem*sizeof(float));
    vel_r_unprc=malloc(elem*sizeof(float));
    vel_theta_unprc=malloc(elem*sizeof(float));
    pres_unprc=malloc(elem*sizeof(float));
    vel_phi_unprc=malloc(elem*sizeof(float));
    
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid, dont need anymore so just save to dummy variable
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(dens_unprc, sizeof(float),elem, hydroPtr); //data
    fclose(hydroPtr);
    
    /*
    for (i=0;i<R_DIM*THETA_DIM*PHI_DIM;i++)
    {
        if ((i>98784000-5) || (i<5))
        {
            fprintf(fPtr,"Density %d: %0.7e\n", i, *(dens_unprc+i));
            fflush(fPtr);
        }
    }
    */
    
    //velocities divided by c
    //v_r
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 2,"-" );
    modifyFlashName(file_num, hydrofile, frame,1);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    /*
    fprintf(fPtr,"Reading v_r: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid, dont need anymore so just save to dummy variable
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(vel_r_unprc, sizeof(float),elem, hydroPtr);
    fclose(hydroPtr);
    
    /*
    for (i=0;i<5;i++)
    {
        fprintf(fPtr,"V_r %d: %e\n", i, *(vel_r_unprc+i));
        fflush(fPtr);
    }
    */
     
    //v_theta
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 3,"-" );
    modifyFlashName(file_num, hydrofile, frame,1);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    /*
    fprintf(fPtr,"Reading v_theta: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid, dont need anymore so just save to dummy variable
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(vel_theta_unprc, sizeof(float),elem, hydroPtr);
    fclose(hydroPtr);
    
    /*
    for (i=0;i<5;i++)
    {
        fprintf(fPtr,"V_theta %d: %e\n", i, *(vel_theta_unprc+i));
        fflush(fPtr);
    }
    */
    
    //v_phi
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 4,"-" );
    modifyFlashName(file_num, hydrofile, frame,1);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    /*
    fprintf(fPtr,"Reading v_phi: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid, dont need anymore so just save to dummy variable
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(vel_phi_unprc, sizeof(float),elem, hydroPtr);
    fclose(hydroPtr);
    
    /*
    for (i=0;i<5;i++)
    {
        fprintf(fPtr,"V_phi %d: %e\n", i, *(vel_phi_unprc+i));
        fflush(fPtr);
    }
    */
    
    //pressure (divided by c^2)
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 8,"-" );
    modifyFlashName(file_num, hydrofile, frame,1);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    /*
    fprintf(fPtr,"Reading pres: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid, dont need anymore so just save to dummy variable
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(pres_unprc, sizeof(float),elem, hydroPtr);
    fclose(hydroPtr);
    /*
    for (i=PHI_DIM-1;i<PHI_DIM;i++)
    {
        for (j=THETA_DIM-1;j<THETA_DIM;j++)
        {
            for (k=R_DIM-5;k<R_DIM;k++)
            {
        
            fprintf(fPtr,"Pres %d: %e\n", (i*R_DIM*THETA_DIM + j*R_DIM + k  ), *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )));
            fflush(fPtr);
        
            }
        }
    }
    */
    
    // see how many elements there are to test if reading correctly
    /*
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    while (1 == fread(&buffer,sizeof(float),1,hydroPtr))
    {
        elem++;
    }
    //fread(pres_unprc, sizeof(double)*HYDRO_DIM,HYDRO_DIM, hydroPtr);
    fclose(hydroPtr);
    
    fprintf(fPtr,"Elem %d\n", elem);
    */
    
    //R
    //remapping_indexes=getIndexesForRadialRemapping(hydro_prefix); //can run this once in debug mode to find out delta index for each remapping and number of total r elements 
    //for given set of remappings on July 12th 2017, grid00-x1.data[420]=grid01-x1.data[0], grid01-x1.data[420]=grid02-x1.data[0], etc. total number of r is 3780
    
    
    if (frame<=1300)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 0,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+0))+r_min_index ; //this if I am not hardcoding the dr index values
        adjusted_remapping_index=(0*420)+r_min_index;
    }
    else if (frame<=2000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 1,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+1))+r_min_index;
        adjusted_remapping_index=(1*420)+r_min_index;
    }
    else if (frame<=10000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 2,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+2))+r_min_index;
        adjusted_remapping_index=(2*420)+r_min_index;
    }
    else if (frame<=20000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 3,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+3))+r_min_index;
        adjusted_remapping_index=(3*420)+r_min_index;
    }
    else if (frame<=35000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 4,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+4))+r_min_index;
        adjusted_remapping_index=(4*420)+r_min_index;
    }
    else if (frame<=50000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 5,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+5))+r_min_index;
        adjusted_remapping_index=(5*420)+r_min_index;
    }
    else if (frame<=60000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 6,"-x1.data" );
        //adjusted_remapping_index=(*(remapping_indexes+6))+r_min_index;
        adjusted_remapping_index=(6*420)+r_min_index;
    }
    
    fprintf(fPtr,"Reading Radius: %s\n", hydrofile);
    fflush(fPtr);
    
    hydroPtr=fopen(hydrofile, "r");
    
    i=0;
    while (i<R_DIM)
    {
        fscanf(hydroPtr, "%lf", (r_unprc+i));  //read value
        fgets(buf, 3,hydroPtr); //read comma
        /*
        if (i<5)
        {
            fprintf(fPtr,"R %d: %e\n", i, *(r_unprc+i));
            fflush(fPtr);
        }
        */
        i++;
    }
    
    fclose(hydroPtr);
    
    r_edge=malloc(sizeof(double)*(3780+1));
    dr=malloc(sizeof(double)*(3780));
    
    //calculate radial grid edges
    *(r_edge+0)=r_in;
    i=0;
    for (i=1;i<3780;i++)
    {
        *(r_edge+i)=(*(r_edge+i-1))+((*(r_edge+i-1))*(M_PI/560)/(1+((*(r_edge+i-1))/r_ref))); //r_i = r_(i-1) + Dq r_(i-1) [1 + r_(i-1)/r0]-1
        *(dr+i-1)=(*(r_edge+i))-(*(r_edge+i-1));
        /*
        if (i<5)
        {
            fprintf(fPtr,"R Edge: %d: %e Dr: %e\n", i, *(r_edge+i), *(dr+i-1));
            fflush(fPtr);
        }
        */
    }
    free(r_edge);
    
    //Theta
    snprintf(hydrofile,sizeof(hydrofile),"%s%s",hydro_prefix,"grid-x2.data" );
    
    fprintf(fPtr,"Reading Theta: %s\n", hydrofile);
    fflush(fPtr);
    
    hydroPtr=fopen(hydrofile, "r");
    
    i=0;
    while (i<THETA_DIM)
    {
        fscanf(hydroPtr, "%lf", (theta_unprc+i));  //read value
        fgets(buf, 3,hydroPtr); //read comma
        /*
        if (i<5)
        {
            fprintf(fPtr,"R %d: %e\n", i, *(theta_unprc+i));
            fflush(fPtr);
        }
        */
        i++;
    }
    
    fclose(hydroPtr);
    
    //Phi
    snprintf(hydrofile,sizeof(hydrofile),"%s%s",hydro_prefix,"grid-x3.data" );
    /*
    fprintf(fPtr,"Reading Phi: %s\n", hydrofile);
    fflush(fPtr);
    */
    hydroPtr=fopen(hydrofile, "r");
    
    i=0;
    while (i<PHI_DIM)
    {
        fscanf(hydroPtr, "%lf", (phi_unprc+i));  //read value
        fgets(buf, 3,hydroPtr); //read comma
        /*
        if (i<5)
        {
            fprintf(fPtr,"R %d: %e\n", i, *(phi_unprc+i));
            fflush(fPtr);
        }
        */
        i++;
    }
    
    fclose(hydroPtr);
    
    //limit number of array elements PUT WHILE LOOP TO MAKE SURE NUMBER OF ELEMENTS >0
    elem_factor=0;
    elem=0;
    while (elem==0)
    {
        elem=0;
        elem_factor++;
        for (i=0;i<(phi_max_index+1-phi_min_index);i++)
        {
            for (j=0;j<(theta_max_index+1-theta_min_index);j++)
            {
                for (k=0;k<(r_max_index+1-r_min_index);k++)
                {
                    r_index=r_min_index+k;
                    //if I have photons do selection differently than if injecting photons
                    if (ph_inj_switch==0)
                    {
                        //printf("R's:%d, %e\n", k, *(r_unprc+r_index));
                        //if calling this function when propagating photons, choose blocks based on where the photons are
                        if (((ph_rmin - elem_factor*C_LIGHT/fps)<(*(r_unprc+r_index))) && (*(r_unprc+r_index)  < (ph_rmax + elem_factor*C_LIGHT/fps) ))
                        {
                            // *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )
                            elem++;
                        }
                    }
                    else
                    {
                        //if calling this function to inject photons choose blocks based on injection parameters, r_inj, which is sufficient 
                        if (((r_inj - elem_factor*C_LIGHT/fps)<(*(r_unprc+r_index))) && (*(r_unprc+r_index)  < (r_inj + elem_factor*C_LIGHT/fps) ))
                        {
                            // *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )
                            elem++;
                        }
                    
                    }
        
                }
            }
        }
    }
    
    fprintf(fPtr,"Number of post restricted Elems: %d %e\n", elem, r_inj);
    //fprintf(fPtr,"Ph_min, Ph_max: %e, %e\n With c: min: %e max: %e \n", ph_rmin, ph_rmax, (ph_rmin - 2*C_LIGHT/fps), (ph_rmax + 2*C_LIGHT/fps));
    fflush(fPtr);
    
    //allocate space for new set of data
    (*pres)=malloc (elem * sizeof (double ));
    (*velx)=malloc (elem * sizeof (double ));
    (*vely)=malloc (elem * sizeof (double ));
    (*velz)=malloc (elem * sizeof (double ));
    (*dens)=malloc (elem * sizeof (double ));
    (*x)=malloc (elem * sizeof (double ));
    (*y)=malloc (elem * sizeof (double ));
    (*z)=malloc (elem * sizeof (double ));
    (*r)=malloc (elem * sizeof (double ));
    (*theta)=malloc (elem * sizeof (double ));
    (*phi)=malloc (elem * sizeof (double ));
    (*gamma)=malloc (elem * sizeof (double ));
    (*dens_lab)=malloc (elem * sizeof (double ));
    (*szx)=malloc (elem * sizeof (double )); //theta and phi resolution 
    (*szy)=malloc (elem * sizeof (double )); //r resolution
    (*temp)=malloc (elem * sizeof (double ));
    
    //limit number of array elements
    elem=0;
    for (i=0;i<(phi_max_index+1-phi_min_index);i++)
    {
        for (j=0;j<(theta_max_index+1-theta_min_index);j++)
        {
            for (k=0;k<(r_max_index+1-r_min_index);k++)
            {
                r_index=r_min_index+k; //look at indexes of r that are included in small hydro file
                theta_index=theta_min_index+j;
                phi_index=phi_min_index+i;
                dr_index=adjusted_remapping_index+k;
                hydro_index=(i*(r_max_index+1-r_min_index)*(theta_max_index+1-theta_min_index) + j*(r_max_index+1-r_min_index) + k  );
                
                //if I have photons do selection differently than if injecting photons
                if (ph_inj_switch==0)
                {
                    //if calling this function when propagating photons, choose blocks based on where the photons are
                    if (((ph_rmin - elem_factor*C_LIGHT/fps)<(*(r_unprc+r_index))) && (*(r_unprc+r_index)  < (ph_rmax + elem_factor*C_LIGHT/fps) ))
                    {
                        (*pres)[elem] = *(pres_unprc+hydro_index);
                        (*dens)[elem] = *(dens_unprc+hydro_index);
                        (*temp)[elem] =  pow(3*(*(pres_unprc+hydro_index))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                        
                        (*gamma)[elem] = pow(pow(1.0-(pow(*(vel_r_unprc+hydro_index),2)+ pow(*(vel_theta_unprc+hydro_index),2)+pow(*(vel_phi_unprc+hydro_index),2)),0.5),-1);
                        (*dens_lab)[elem] = (*(dens_unprc+hydro_index))*pow(pow(1.0-(pow(*(vel_r_unprc+hydro_index),2)+ pow(*(vel_theta_unprc+hydro_index),2)+pow(*(vel_phi_unprc+hydro_index),2)),0.5),-1);
                        (*r)[elem] = *(r_unprc+r_index);
                        (*theta)[elem] = *(theta_unprc+theta_index);
                        (*phi)[elem] = *(phi_unprc+phi_index);
                        (*x)[elem] = (*(r_unprc+r_index))*sin(*(theta_unprc+theta_index))*cos(*(phi_unprc+phi_index));
                        (*y)[elem] = (*(r_unprc+r_index))*sin(*(theta_unprc+theta_index))*sin(*(phi_unprc+phi_index));
                        (*z)[elem] = (*(r_unprc+r_index))*cos(*(theta_unprc+theta_index));
                        (*szx)[elem] =  *(dr+dr_index);
                        (*szy)[elem] =  M_PI/560;
                        (*velx)[elem]=((*(vel_r_unprc+hydro_index))*sin(*(theta_unprc+theta_index))*cos(*(phi_unprc+phi_index))) + ((*(vel_theta_unprc+hydro_index))*cos(*(theta_unprc+theta_index))*cos(*(phi_unprc+phi_index))) - ((*(vel_phi_unprc+hydro_index))*sin(*(phi_unprc+phi_index)));
                        (*vely)[elem]=((*(vel_r_unprc+hydro_index))*sin(*(theta_unprc+theta_index))*sin(*(phi_unprc+phi_index))) + ((*(vel_theta_unprc+hydro_index))*cos(*(theta_unprc+theta_index))*sin(*(phi_unprc+phi_index))) + ((*(vel_phi_unprc+hydro_index))*cos(*(phi_unprc+phi_index)));
                        (*velz)[elem]=((*(vel_r_unprc+hydro_index))*cos(*(theta_unprc+theta_index))) - ((*(vel_theta_unprc+hydro_index))*sin(*(theta_unprc+theta_index)));

                        elem++;
                        
                    }
                }
                else
                {
                    //if calling this function to inject photons choose blocks based on injection parameters, r_inj, which is sufficient 
                    if (((r_inj - elem_factor*C_LIGHT/fps)<(*(r_unprc+r_index))) && (*(r_unprc+r_index)  < (r_inj + elem_factor*C_LIGHT/fps) ))
                    {
                        (*pres)[elem] = *(pres_unprc+hydro_index);
                        (*dens)[elem] = *(dens_unprc+hydro_index);
                        (*temp)[elem] =  pow(3*(*(pres_unprc+hydro_index))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                        (*gamma)[elem] = pow(pow(1.0-(pow(*(vel_r_unprc+hydro_index),2)+ pow(*(vel_theta_unprc+hydro_index),2)+pow(*(vel_phi_unprc+hydro_index),2)),0.5),-1);
                        (*dens_lab)[elem] = (*(dens_unprc+hydro_index))*pow(pow(1.0-(pow(*(vel_r_unprc+hydro_index),2)+ pow(*(vel_theta_unprc+hydro_index),2)+pow(*(vel_phi_unprc+hydro_index),2)),0.5),-1);
                        (*r)[elem] = *(r_unprc+r_index);
                        (*theta)[elem] = *(theta_unprc+theta_index);
                        (*phi)[elem] = *(phi_unprc+phi_index);
                        (*x)[elem] = (*(r_unprc+r_index))*sin(*(theta_unprc+theta_index))*cos(*(phi_unprc+phi_index));
                        (*y)[elem] = (*(r_unprc+r_index))*sin(*(theta_unprc+theta_index))*sin(*(phi_unprc+phi_index));
                        (*z)[elem] = (*(r_unprc+r_index))*cos(*(theta_unprc+theta_index));
                        (*szx)[elem] =  *(dr+dr_index);
                        (*szy)[elem] =  M_PI/560;
                        (*velx)[elem]=((*(vel_r_unprc+hydro_index))*sin(*(theta_unprc+theta_index))*cos(*(phi_unprc+phi_index))) + ((*(vel_theta_unprc+hydro_index))*cos(*(theta_unprc+theta_index))*cos(*(phi_unprc+phi_index))) - ((*(vel_phi_unprc+hydro_index))*sin(*(phi_unprc+phi_index)));
                        (*vely)[elem]=((*(vel_r_unprc+hydro_index))*sin(*(theta_unprc+theta_index))*sin(*(phi_unprc+phi_index))) + ((*(vel_theta_unprc+hydro_index))*cos(*(theta_unprc+theta_index))*sin(*(phi_unprc+phi_index))) + ((*(vel_phi_unprc+hydro_index))*cos(*(phi_unprc+phi_index)));
                        (*velz)[elem]=((*(vel_r_unprc+hydro_index))*cos(*(theta_unprc+theta_index))) - ((*(vel_theta_unprc+hydro_index))*sin(*(theta_unprc+theta_index)));

                        elem++;
                    }
                    
                }
        
            }
        }
    }
    
    *number=elem;
    
    free(pres_unprc); free(dens_unprc); free(r_unprc); free(theta_unprc); free(phi_unprc);free(dr);free(vel_r_unprc); free(vel_theta_unprc); free(vel_phi_unprc); 
    
}


void photonInjection3D( struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, int array_length, double fps, double theta_min, double theta_max,\
 double *x, double *y, double *z, double *szx, double *szy, double *r, double *theta, double *phi, double *temps, double *vx, double *vy, double *vz, gsl_rng * rand, FILE *fPtr)
{
    int i=0, block_cnt=0, *ph_dens=NULL, ph_tot=0, j=0,k=0;
    double ph_dens_calc=0.0, fr_dum=0.0, y_dum=0.0, yfr_dum=0.0, fr_max=0, bb_norm=0, position_phi, ph_weight_adjusted, theta_prime=0;
    double com_v_phi, com_v_theta, *p_comv=NULL, *boost=NULL; //comoving phi, theta, comoving 4 momentum for a photon, and boost for photon(to go to lab frame)
    double *l_boost=NULL; //pointer to hold array of lorentz boost, to lab frame, values
    float num_dens_coeff;
    
    if (spect=='w') //from MCRAT paper, w for wien spectrum 
    {
        num_dens_coeff=8.44;
        //printf("in wien spectrum\n");
    }
    else
    {
        num_dens_coeff=20.29; //this is for black body spectrum
        //printf("in BB spectrum");
    }
    
    //find how many blocks are near the injection radius within the angles defined in mc.par, get temperatures and calculate number of photons to allocate memory for 
    //and then rcord which blocks have to have "x" amount of photons injected there
    printf("%e, %e\n",*(phi+i), theta_max);
    for(i=0;i<array_length;i++)
    {
        //look at all boxes in width delta r=c/fps and within angles we are interested in NEED TO modify for RIKEN data- dont need r anymore, just theta and phi? (didnt work), just look at pojection on x-z plane
        theta_prime=acos(*(y+i)/(*(r+i))); //jet axis here is the y axis
            if ( (theta_prime< theta_max) && (theta_prime >= theta_min) ) //(*(r+i) > (r_inj - C_LIGHT/fps))  &&   (*(r+i)  < (r_inj + C_LIGHT/fps)  ) &&
            {
                //printf("%e\n", theta_prime );
                block_cnt++;
            }
    }
    printf("Blocks: %d\n", block_cnt);
    
    ph_dens=malloc(block_cnt * sizeof(int));
    
    //calculate the photon density for each block and save it to the array
    j=0;
    ph_tot=0;
    ph_weight_adjusted=ph_weight;
    //printf("%d %d\n", max_photons, min_photons);
    while ((ph_tot>max_photons) || (ph_tot<min_photons) )
    {
        j=0;
        ph_tot=0;
        //allocate memory to record density of photons for each block
        //ph_dens=malloc(block_cnt * sizeof(int));
        
        for (i=0;i<array_length;i++)
        {
            //printf("%d\n",i);
            //printf("%e, %e, %e, %e, %e, %e\n", *(r+i),(r_inj - C_LIGHT/fps), (r_inj + C_LIGHT/fps), *(theta+i) , theta_max, theta_min);
                //NEED TO modify for RIKEN data - modified
                theta_prime=acos(*(y+i)/(*(r+i)));
                if ( (theta_prime< theta_max) && (theta_prime >= theta_min) )
                {
                    //NEED TO modify for RIKEN data - modified
                    ph_dens_calc=(num_dens_coeff*pow(*(temps+i),3.0)*pow(*(r+i),2)*sin(*(theta+i))* pow(*(szy+i),2.0)*(*(szx+i)) /(ph_weight_adjusted))*pow(pow(1.0-(pow(*(vx+i),2)+pow(*(vy+i),2)+pow(*(vz+i),2)),0.5),-1) ; //a*T^3/(weight) dV, dV=2*PI*x*dx^2,
                     
                     (*(ph_dens+j))=gsl_ran_poisson(rand,ph_dens_calc) ; //choose from poission distribution with mean of ph_dens_calc
                     
                    //printf("%d, %lf \n",*(ph_dens+j), ph_dens_calc);
                    
                     //sum up all the densities to get total number of photons
                     ph_tot+=(*(ph_dens+j));
                     
                     j++;
                }
        }

        if (ph_tot>max_photons)
        {
            //if the number of photons is too big make ph_weight larger
            ph_weight_adjusted*=10;
            //free(ph_dens);
        }
        else if (ph_tot<min_photons)
        {
            ph_weight_adjusted*=0.5;
            //free(ph_dens);
        }
        //printf("dens: %d, photons: %d\n", *(ph_dens+(j-1)), ph_tot);
         
    }
        
    printf("%d\n", ph_tot);
    
    //allocate memory for that many photons and also allocate memory to hold comoving 4 momentum of each photon and the velocity of the fluid
    (*ph)=malloc (ph_tot * sizeof (struct photon ));
    
    p_comv=malloc(4*sizeof(double));
    boost=malloc(3*sizeof(double));
    l_boost=malloc(4*sizeof(double));
    
    
    //go through blocks and assign random energies/locations to proper number of photons
    ph_tot=0;
    k=0;
    for (i=0;i<array_length;i++)
    {
        theta_prime=acos(*(y+i)/(*(r+i)));
        if ( (theta_prime< theta_max) && (theta_prime >= theta_min) )  //NEED TO modify for RIKEN data - modified
        {

            //*(temps+i)=0.76*(*(temps+i));
            for(j=0;j<( *(ph_dens+k) ); j++ )
            {
                    //have to get random frequency for the photon comoving frequency
                    y_dum=1; //initalize loop
                    yfr_dum=0;
                    while (y_dum>yfr_dum)
                    {
                        fr_dum=gsl_rng_uniform_pos(rand)*6.3e11*(*(temps+i)); //in Hz
                        //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
                        y_dum=gsl_rng_uniform_pos(rand);
                        //printf("%lf ",fr_dum);
                        
                        if (spect=='w')
                        {
                            yfr_dum=(1.0/(1.29e31))*pow((fr_dum/(*(temps+i))),3.0)/(exp((PL_CONST*fr_dum)/(K_B*(*(temps+i)) ))-1); //curve is normalized to maximum
                        }
                        else
                        {
                            fr_max=(5.88e10)*(*(temps+i));//(C_LIGHT*(*(temps+i)))/(0.29); //max frequency of bb
                            bb_norm=(PL_CONST*fr_max * pow((fr_max/C_LIGHT),2.0))/(exp(PL_CONST*fr_max/(K_B*(*(temps+i))))-1); //find value of bb at fr_max
                            yfr_dum=((1.0/bb_norm)*PL_CONST*fr_dum * pow((fr_dum/C_LIGHT),2.0))/(exp(PL_CONST*fr_dum/(K_B*(*(temps+i))))-1); //curve is normalized to vaue of bb @ max frequency
                        }
                        //printf("%lf, %lf,%lf,%e \n",(*(temps+i)),fr_dum, y_dum, yfr_dum);
                        
                    }
                   //printf("%lf\n ",fr_dum);
                   //position_phi= gsl_rng_uniform(rand)*2*M_PI; //NEED TO modify for RIKEN data-modified, dont need anymore
                   com_v_phi=gsl_rng_uniform(rand)*2*M_PI; 
                   com_v_theta=acos((gsl_rng_uniform(rand)*2)-1);
                   //printf("%lf, %lf, %lf\n", position_phi, com_v_phi, com_v_theta);
                   
                   //populate 4 momentum comoving array
                   *(p_comv+0)=PL_CONST*fr_dum/C_LIGHT;
                   *(p_comv+1)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*cos(com_v_phi);
                   *(p_comv+2)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*sin(com_v_phi);
                   *(p_comv+3)=(PL_CONST*fr_dum/C_LIGHT)*cos(com_v_theta);
                   
                    //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
                    //NEED TO modify for RIKEN data - modified
                    *(boost+0)=-1*(*(vx+i));
                    *(boost+1)=-1*(*(vy+i));
                    *(boost+2)=-1*(*(vz+i));
                    //printf("%lf, %lf, %lf\n", *(boost+0), *(boost+1), *(boost+2));
                    
                    //boost to lab frame
                    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
                    //printf("Assignemnt: %e, %e, %e, %e\n", *(l_boost+0), *(l_boost+1), *(l_boost+2),*(l_boost+3));
                   
                (*ph)[ph_tot].p0=(*(l_boost+0));
                (*ph)[ph_tot].p1=(*(l_boost+1));
                (*ph)[ph_tot].p2=(*(l_boost+2));
                (*ph)[ph_tot].p3=(*(l_boost+3));
                //NEED TO modify for RIKEN data-modified
                (*ph)[ph_tot].r0= (*(x+i)); //put photons @ center of box that they are supposed to be in with random phi 
                (*ph)[ph_tot].r1=(*(y+i)) ;
                (*ph)[ph_tot].r2=(*(z+i)); //y coordinate in flash becomes z coordinate in MCRaT
                (*ph)[ph_tot].num_scatt=0;
                (*ph)[ph_tot].weight=ph_weight_adjusted;
                //printf("%d\n",ph_tot);
                ph_tot++;
            }
            k++;
        }
    }
    
    *ph_num=ph_tot; //save number of photons
    //printf(" %d: %d\n", *(ph_dens+(k-1)), *ph_num);
    free(ph_dens); free(p_comv);free(boost); free(l_boost);
    
}

void phMinMax(struct photon *ph, int ph_num, double *min, double *max, double *min_theta, double *max_theta, FILE *fPtr)
{
    double temp_r_max=0, temp_r_min=DBL_MAX, temp_theta_max=0, temp_theta_min=DBL_MAX;
    int i=0, num_thread=omp_get_num_threads();
    double ph_r=0, ph_theta=0;
    
#pragma omp parallel for num_threads(num_thread) firstprivate(ph_r, ph_theta) reduction(min:temp_r_min) reduction(max:temp_r_max) reduction(min:temp_theta_min) reduction(max:temp_theta_max)
    for (i=0;i<ph_num;i++)
    {
        if ((ph+i)->weight != 0)
        {
            ph_r=pow(pow( ((ph+i)->r0), 2.0) + pow(((ph+i)->r1),2.0 ) + pow(((ph+i)->r2) , 2.0),0.5);
            ph_theta=acos(((ph+i)->r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            if (ph_r > temp_r_max )
            {
                temp_r_max=ph_r;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_r<temp_r_min)
            {
                temp_r_min=ph_r;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            if (ph_theta > temp_theta_max )
            {
                temp_theta_max=ph_theta;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_theta<temp_theta_min)
            {
                temp_theta_min=ph_theta;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
        }
    }
    
    *max=temp_r_max;
    *min=temp_r_min;
    *max_theta=temp_theta_max;
    *min_theta=temp_theta_min;
}

int *getIndexesForRadialRemapping(char hydro_prefix[200])
{
    FILE *hydroPtr=NULL;
    char hydrofile[200]="";
    char buf[10]="";
    int i=0, j=0;
    int *remapping_start_index=malloc(sizeof(int)*7); //what index out of total range of r does each remapping begin at
    double r_in=1e10, r_ref=2e13;
    double *r_unprc_0=malloc(sizeof(double)*R_DIM), *r_unprc_1=malloc(sizeof(double)*R_DIM), *r_unprc_2=malloc(sizeof(double)*R_DIM), *r_unprc_3=malloc(sizeof(double)*R_DIM);
    double *r_unprc_4=malloc(sizeof(double)*R_DIM), *r_unprc_5=malloc(sizeof(double)*R_DIM), *r_unprc_6=malloc(sizeof(double)*R_DIM);
    double *r_edge=NULL, *dr=NULL, *rPtr=NULL;
    
    for (i=0;i<7;i++)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", i,"-x1.data" );
    
        hydroPtr=fopen(hydrofile, "r");
    
        j=0;
        while (j<R_DIM)
        {
            switch (i)
            {
                case 0: fscanf(hydroPtr, "%lf", (r_unprc_0+j));  //read value
                case 1: fscanf(hydroPtr, "%lf", (r_unprc_1+j));
                case 2: fscanf(hydroPtr, "%lf", (r_unprc_2+j));
                case 3: fscanf(hydroPtr, "%lf", (r_unprc_3+j));
                case 4: fscanf(hydroPtr, "%lf", (r_unprc_4+j));
                case 5: fscanf(hydroPtr, "%lf", (r_unprc_5+j));
                case 6: fscanf(hydroPtr, "%lf", (r_unprc_6+j));
            
            }
            fgets(buf, 3,hydroPtr); //read comma
            /*
            if (i<5)
            {
                fprintf(fPtr,"R %d: %e\n", i, *(r_unprc+i));
                fflush(fPtr);
            }
            */
            j++;
        }
    
        fclose(hydroPtr);
    }
    
    //calculate the indexes in which each remapping takes over
    j=0; //keeps track of indexes of all of the possible r values
    i=0; //keeps track of index within a certain remapping, when get to R_DIM know that were @ end of last remapping b/c remappings overlap with one another
    rPtr=r_unprc_0; //start off looking at 0th remapping
    *(remapping_start_index+1)=j; //0th remapping starts at index 0
    while (i<R_DIM)
    {
        if (*(rPtr+i)== *(r_unprc_1+0))
        {
            //if the element of the 0th remapping is equal to the 1st element of the 1st remapping, start to look at the 1st remapping
            rPtr=r_unprc_1;
            i=0;
            *(remapping_start_index+1)=j;
        }
        else if (*(rPtr+i)== *(r_unprc_2+0))
        {
            rPtr=r_unprc_2;
            i=0;
            *(remapping_start_index+2)=j;
        }
        else if (*(rPtr+i)== *(r_unprc_3+0))
        {
            rPtr=r_unprc_3;
            i=0;
            *(remapping_start_index+3)=j;
        }
        else if (*(rPtr+i)== *(r_unprc_4+0))
        {
            rPtr=r_unprc_4;
            i=0;
            *(remapping_start_index+4)=j;
        }
        else if (*(rPtr+i)== *(r_unprc_5+0))
        {
            rPtr=r_unprc_5;
            i=0;
            *(remapping_start_index+5)=j;
        }
        else if (*(rPtr+i)== *(r_unprc_6+0))
        {
            rPtr=r_unprc_6;
            i=0;
            *(remapping_start_index+6)=j;
        }
        
        j++;
        i++;
    }
    
    printf("Indexes %d, %d, %d, %d, %d, %d, %d\n Elems: %d\n", *(remapping_start_index+0), *(remapping_start_index+1), *(remapping_start_index+2), *(remapping_start_index+3), *(remapping_start_index+4), *(remapping_start_index+5), *(remapping_start_index+6), j);
    //exit(0);
    r_edge=malloc(sizeof(double)*(j+1));
    dr=malloc(sizeof(double)*j);
    
    //calculate radial grid edges
    *(r_edge+0)=r_in;
    i=0;
    for (i=1;i<j;i++)
    {
        *(r_edge+i)=(*(r_edge+i-1))+((*(r_edge+i-1))*(M_PI/560)/(1+((*(r_edge+i-1))/r_ref))); //r_i = r_(i-1) + Dq r_(i-1) [1 + r_(i-1)/r0]-1
        *(dr+i-1)=(*(r_edge+i))-(*(r_edge+i-1));
        /*
        if (i<5)
        {
            fprintf(fPtr,"R Edge: %d: %e Dr: %e\n", i, *(r_edge+i), *(dr+i-1));
            fflush(fPtr);
        }
        */
    }
    free(r_edge);
    free(r_unprc_0);
    free(r_unprc_1);
    free(r_unprc_2);
    free(r_unprc_3);
    free(r_unprc_4);
    free(r_unprc_5);
    free(r_unprc_6);
    
    return remapping_start_index;
    
}
