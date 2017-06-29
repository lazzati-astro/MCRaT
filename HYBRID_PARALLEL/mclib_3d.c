#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
    char hydrofile[200]="", file_num[200]="", full_file[200]=""  ;
    char buf[10]="";
    int i=0, j=0, k=0, elem=0;
    float buffer=0;
    float *dens_unprc=malloc(sizeof(float)*R_DIM*THETA_DIM*PHI_DIM);
    float *vel_r_unprc=malloc(sizeof(float)*R_DIM*THETA_DIM*PHI_DIM);
    float *vel_theta_unprc=malloc(sizeof(float)*R_DIM*THETA_DIM*PHI_DIM);
    float *vel_phi_unprc=malloc(sizeof(float)*R_DIM*THETA_DIM*PHI_DIM);
    float *pres_unprc=malloc(sizeof(float)*R_DIM*THETA_DIM*PHI_DIM);
    double ph_rmin=0, ph_rmax=0;
    double r_in=1e10, r_ref=2e13;
    double *r_edge=malloc(sizeof(double)*(R_DIM+1));
    double *dr=malloc(sizeof(double)*(R_DIM));
    double *r_unprc=malloc(sizeof(double)*R_DIM);
    double *theta_unprc=malloc(sizeof(double)*THETA_DIM);
    double *phi_unprc=malloc(sizeof(double)*PHI_DIM);
    
    if (ph_inj_switch==0)
    {
        ph_rmin=min_r;
        ph_rmax=max_r;
    }
    
    
    
    //density
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 1,"-" );
    modifyFlashName(file_num, hydrofile, frame,1);
    
    fprintf(fPtr,">> Opening file %s\n", file_num);
    fflush(fPtr);
    
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, ".data");
    /*
    fprintf(fPtr,"Reading Density: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(dens_unprc, sizeof(float)*R_DIM*THETA_DIM*PHI_DIM,R_DIM*THETA_DIM*PHI_DIM, hydroPtr); //data
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
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, ".data");
    /*
    fprintf(fPtr,"Reading v_r: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(vel_r_unprc, sizeof(float)*R_DIM*THETA_DIM*PHI_DIM,R_DIM*THETA_DIM*PHI_DIM, hydroPtr);
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
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, ".data");
    /*
    fprintf(fPtr,"Reading v_theta: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(vel_theta_unprc, sizeof(float)*R_DIM*THETA_DIM*PHI_DIM,R_DIM*THETA_DIM*PHI_DIM, hydroPtr);
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
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, ".data");
    /*
    fprintf(fPtr,"Reading v_phi: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(vel_phi_unprc, sizeof(float)*R_DIM*THETA_DIM*PHI_DIM,R_DIM*THETA_DIM*PHI_DIM, hydroPtr);
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
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, ".data");
    /*
    fprintf(fPtr,"Reading pres: %s\n", full_file);
    fflush(fPtr);
    */
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran 
    fread(pres_unprc, sizeof(float)*R_DIM*THETA_DIM*PHI_DIM,R_DIM*THETA_DIM*PHI_DIM, hydroPtr);
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
    if (frame<=1300)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 0,"-x1.data" );
    }
    else if (frame<=2000)
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 1,"-x1.data" );
    }
    else
    {
        snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"grid0", 2,"-x1.data" );
    }
    /*
    fprintf(fPtr,"Reading Radius: %s\n", hydrofile);
    fflush(fPtr);
    */
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
    
    //calculate radial grid edges
    *(r_edge+0)=r_in;
    i=0;
    for (i=1;i<R_DIM;i++)
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
    /*
    fprintf(fPtr,"Reading Theta: %s\n", hydrofile);
    fflush(fPtr);
    */
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
    
    //limit number of array elements
    elem=0;
    for (i=0;i<PHI_DIM;i++)
    {
        for (j=0;j<THETA_DIM;j++)
        {
            for (k=0;k<R_DIM;k++)
            {
                //if I have photons do selection differently than if injecting photons
                if (ph_inj_switch==0)
                {
                    //if calling this function when propagating photons, choose blocks based on where the photons are
                    if (((ph_rmin - 2*C_LIGHT/fps)<(*(r_unprc+k))) && (*(r_unprc+k)  < (ph_rmax + 2*C_LIGHT/fps) ))
                    {
                        // *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )
                        elem++;
                    }
                }
                else
                {
                    //if calling this function to inject photons choose blocks based on injection parameters, r_inj, which is sufficient 
                    if (((r_inj - C_LIGHT/fps)<(*(r_unprc+k))) && (*(r_unprc+k)  < (r_inj + C_LIGHT/fps) ))
                    {
                        // *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )
                        elem++;
                    }
                    
                }
        
            }
        }
    }
    
    //fprintf(fPtr,"Number of post restricted Elems: %d %e\n", elem, r_inj);
    
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
    for (i=0;i<PHI_DIM;i++)
    {
        for (j=0;j<THETA_DIM;j++)
        {
            for (k=0;k<R_DIM;k++)
            {
                //if I have photons do selection differently than if injecting photons
                if (ph_inj_switch==0)
                {
                    //if calling this function when propagating photons, choose blocks based on where the photons are
                    if (((ph_rmin - 2*C_LIGHT/fps)<(*(r_unprc+k))) && (*(r_unprc+k)  < (ph_rmax + 2*C_LIGHT/fps) ))
                    {
                        (*pres)[elem] = *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  ));
                        (*dens)[elem] = *(dens_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  ));
                        (*temp)[elem] =  pow(3*(*(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                        
                        (*gamma)[elem] = pow(pow(1.0-(pow(*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+ pow(*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+pow(*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)),0.5),-1);
                        (*dens_lab)[elem] = (*(dens_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )))*pow(pow(1.0-(pow(*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+ pow(*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+pow(*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)),0.5),-1);
                        (*r)[elem] = *(r_unprc+k);
                        (*theta)[elem] = *(theta_unprc+j);
                        (*phi)[elem] = *(phi_unprc+i);
                        (*x)[elem] = (*(r_unprc+k))*sin(*(theta_unprc+j))*cos(*(phi_unprc+i));
                        (*y)[elem] = (*(r_unprc+k))*sin(*(theta_unprc+j))*sin(*(phi_unprc+i));
                        (*z)[elem] = (*(r_unprc+k))*cos(*(theta_unprc+j));
                        (*szx)[elem] = M_PI/560; 
                        (*szy)[elem] =  *(dr+k);
                        (*velx)[elem]=((*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(theta_unprc+j))*cos(*(phi_unprc+i))) + ((*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(theta_unprc+j))*cos(*(phi_unprc+i))) - ((*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(phi_unprc+i)));
                        (*vely)[elem]=((*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(theta_unprc+j))*sin(*(phi_unprc+i))) + ((*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(theta_unprc+j))*sin(*(phi_unprc+i))) + ((*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(phi_unprc+i)));
                        (*velz)[elem]=((*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(theta_unprc+j))) - ((*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(theta_unprc+j)));

                        elem++;
                        
                    }
                }
                else
                {
                    //if calling this function to inject photons choose blocks based on injection parameters, r_inj, which is sufficient 
                    if (((r_inj - C_LIGHT/fps)<(*(r_unprc+k))) && (*(r_unprc+k)  < (r_inj + C_LIGHT/fps) ))
                    {
                        (*pres)[elem] = *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  ));
                        (*dens)[elem] = *(dens_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  ));
                        (*temp)[elem] =  pow(3*(*(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                        (*gamma)[elem] = pow(pow(1.0-(pow(*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+ pow(*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+pow(*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)),0.5),-1);
                        (*dens_lab)[elem] = (*(dens_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )))*pow(pow(1.0-(pow(*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+ pow(*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)+pow(*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )),2)),0.5),-1);
                        (*r)[elem] = *(r_unprc+k);
                        (*theta)[elem] = *(theta_unprc+j);
                        (*phi)[elem] = *(phi_unprc+i);
                        (*x)[elem] = (*(r_unprc+k))*sin(*(theta_unprc+j))*cos(*(phi_unprc+i));
                        (*y)[elem] = (*(r_unprc+k))*sin(*(theta_unprc+j))*sin(*(phi_unprc+i));
                        (*z)[elem] = (*(r_unprc+k))*cos(*(theta_unprc+j));
                        (*szx)[elem] = M_PI/560; 
                        (*szy)[elem] =  *(dr+k);
                        (*velx)[elem]=((*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(theta_unprc+j))*cos(*(phi_unprc+i))) + ((*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(theta_unprc+j))*cos(*(phi_unprc+i))) - ((*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(phi_unprc+i)));
                        (*vely)[elem]=((*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(theta_unprc+j))*sin(*(phi_unprc+i))) + ((*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(theta_unprc+j))*sin(*(phi_unprc+i))) + ((*(vel_phi_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(phi_unprc+i)));
                        (*velz)[elem]=((*(vel_r_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*cos(*(theta_unprc+j))) - ((*(vel_theta_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k )))*sin(*(theta_unprc+j)));

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
            if ( (theta_prime< theta_max) && (theta_prime > theta_min) ) //(*(r+i) > (r_inj - C_LIGHT/fps))  &&   (*(r+i)  < (r_inj + C_LIGHT/fps)  ) &&
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
                if ( (theta_prime< theta_max) && (theta_prime > theta_min) )
                {
                    //NEED TO modify for RIKEN data - modified
                    ph_dens_calc=(num_dens_coeff*pow(*(temps+i),3.0)*pow(*(r+i),2)*sin(*(theta+i))* pow(*(szx+i),2.0)*(*(szy+i)) /(ph_weight_adjusted))*pow(pow(1.0-(pow(*(vx+i),2)+pow(*(vy+i),2)),0.5),-1) ; //a*T^3/(weight) dV, dV=2*PI*x*dx^2,
                     
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
        if ( (theta_prime< theta_max) && (theta_prime > theta_min) )  //NEED TO modify for RIKEN data - modified
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

 void phMinMax(struct photon *ph, int ph_num, double *min, double *max)
 {
      double temp_r_max=0, temp_r_min=-1;
        int i=0;
      double ph_r=0;
      
      for (i=0;i<ph_num;i++)
    {        
        ph_r=pow(pow( ((ph+i)->r0), 2.0) + pow(((ph+i)->r1),2.0 ) + pow(((ph+i)->r2) , 2.0),0.5);
        if (ph_r > temp_r_max )
        {
            temp_r_max=ph_r;
            //printf("The new max is: %e\n", temp_r_max);
        }
        
        if ((i==0) || (ph_r<temp_r_min))
        {
            temp_r_min=ph_r;
            //printf("The new min is: %e\n", temp_r_min);
        }
    }
    
    *max=temp_r_max;
    *min=temp_r_min;
      
 }
