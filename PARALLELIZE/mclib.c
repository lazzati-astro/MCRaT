#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
#include <omp.h>

#define PROP_DIM1 1
#define PROP_DIM2 8
#define PROP_DIM3 8
#define COORD_DIM1 2

//define constants
const double A_RAD=7.56e-15, C_LIGHT=2.99792458e10, PL_CONST=6.6260755e-27;
const double K_B=1.380658e-16, M_P=1.6726231e-24, THOMP_X_SECT=6.65246e-25, M_EL=(1.6726231e-24)/1836 ;

void printPhotons(struct photon *ph, int num_ph, int frame,char dir[200] )
{
    //function to save the photons' positions and 4 momentum
    int i=0;
    char mc_file_p0[200], mc_file_p1[200],mc_file_p2[200], mc_file_p3[200];
    char mc_file_r0[200], mc_file_r1[200], mc_file_r2[200], mc_file_ns[200];
    FILE *fPtr=NULL, *fPtr1=NULL,*fPtr2=NULL,*fPtr3=NULL,*fPtr4=NULL,*fPtr5=NULL,*fPtr6=NULL,*fPtr7=NULL;
    
    //make strings for proper files
    snprintf(mc_file_p0,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_P0.dat" );
    snprintf(mc_file_p1,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_P1.dat" );
    snprintf(mc_file_p2,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_P2.dat" );
    snprintf(mc_file_p3,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_P3.dat" );
    snprintf(mc_file_r0,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_R0.dat" );
    snprintf(mc_file_r1,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_R1.dat" );
    snprintf(mc_file_r2,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_R2.dat" );
    snprintf(mc_file_ns,sizeof(mc_file_p0),"%s%s%d%s",dir,"mcdata_", frame,"_NS.dat" ); //for number of scatterings each photon went through
    
    //save the energy
    fPtr=fopen(mc_file_p0, "a");
    fPtr1=fopen(mc_file_p1, "a");
    fPtr2=fopen(mc_file_p2, "a");
    fPtr3=fopen(mc_file_p3, "a");
    fPtr4=fopen(mc_file_r0, "a");
    fPtr5=fopen(mc_file_r1, "a");
    fPtr6=fopen(mc_file_r2, "a");
    fPtr7=fopen(mc_file_ns, "a");
    
    //printf("Writing P0\n");
    for (i=0;i<num_ph;i++)
    {
        fprintf(fPtr,"%0.13e\t",  (ph+i)->p0);
        //printf("%d: %0.13e \n", i, (ph+i)->p0);
        
        fprintf(fPtr1,"%0.13e\t",  (ph+i)->p1);
        //printf("%d: %0.13e \n", i, (ph+i)->p1);
        
        fprintf(fPtr2,"%0.13e\t",  (ph+i)->p2);
        //printf("%d: %0.13e \n", i, (ph+i)->p2);
        
        fprintf(fPtr3,"%0.13e\t",  (ph+i)->p3);
        //printf("%d: %0.13e \n", i, (ph+i)->p3);
        
        fprintf(fPtr4,"%0.13e\t",  (ph+i)->r0);
        //printf("%d: %0.13e \n", i, (ph+i)->r0);
        
        fprintf(fPtr5,"%0.13e\t",  (ph+i)->r1);
        //printf("%d: %0.13e \n", i, (ph+i)->r1);
        
        fprintf(fPtr6,"%0.13e\t",  (ph+i)->r2);
        //printf("%d: %0.13e \n", i, (ph+i)->r2);
        
        //fprintf(fPtr7,"%0.13e\t",  *(ph_num_scatt+i));
        fprintf(fPtr7,"%e\t",  (ph+i)->num_scatt);
        //printf("%d: %0.13e \n", i, (ph+i)->num_scatt);
        
    }    
    fclose(fPtr);
    fclose(fPtr1);
    fclose(fPtr2);
    fclose(fPtr3);
    fclose(fPtr4);
    fclose(fPtr5);
    fclose(fPtr6);
    fclose(fPtr7);
    
    //printf("%s\n%s\n%s\n", mc_file_p0, mc_file_r0, mc_file_ns);
}

void saveCheckpoint(char dir[200], int frame, int scatt_frame, int ph_num,double time_now, struct photon *ph, int last_frame )
{
    //function to save data necessary to restart simulation if it ends
    //need to save all photon data 
    FILE *fPtr=NULL;
    char checkptfile[200];
    char restart;
    int i=0;
    
    snprintf(checkptfile,sizeof(checkptfile),"%s%s",dir,"mc_chkpt.dat" );

    
    fPtr=fopen(checkptfile, "wb");
    
    if (scatt_frame!=last_frame)
    {
        restart='c';
        fwrite(&restart, sizeof(char), 1, fPtr);
        fwrite(&frame, sizeof(int), 1, fPtr);
        fwrite(&scatt_frame, sizeof(int), 1, fPtr);
        fwrite(&time_now, sizeof(double), 1, fPtr);
        fwrite(&ph_num, sizeof(int), 1, fPtr);
        
        //for(i=0;i<ph_num;i++)
        //{
            fwrite((ph), sizeof (struct photon )*ph_num, ph_num, fPtr);
        //}
        
    }
    else
    {
        //just finished last iteration of scatt_frame
        restart='r';
        fwrite(&restart, sizeof(char), 1, fPtr);
        fwrite(&frame, sizeof(int), 1, fPtr);
    }
    fclose(fPtr);
    
}

void readCheckpoint(char dir[200], struct photon **ph, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time )
{
    //function to read in data from checkpoint file
    FILE *fPtr=NULL;
    char checkptfile[200];
    int i=0;
    //int frame, scatt_frame, ph_num, i=0;
    struct photon *phHolder=NULL; //pointer to struct to hold data read in from checkpoint file
    
    snprintf(checkptfile,sizeof(checkptfile),"%s%s",dir,"mc_chkpt.dat" );
    
    fPtr=fopen(checkptfile, "rb");
    
    fread(restart, sizeof(char), 1, fPtr);
    printf("%c\n", *restart);
    fread(framestart, sizeof(int), 1, fPtr);
    printf("%d\n", *framestart);
    
    if((*restart)=='c')
    {
        fread(scatt_framestart, sizeof(int), 1, fPtr);
        *scatt_framestart+=1; //add one to start at the next frame after the siomulation was interrrupted
        printf("%d\n", *scatt_framestart);
        fread(time, sizeof(double), 1, fPtr);
        printf("%e\n", *time);
        fread(ph_num, sizeof(int), 1, fPtr);
        printf("%d\n", *ph_num);
        
        phHolder=malloc(sizeof(struct photon));
        (*ph)=malloc(sizeof(struct photon)*(*ph_num)); //allocate memory to hold photon data
        
        for (i=0;i<(*ph_num);i++)
        {
            fread(phHolder, sizeof(struct photon), 1, fPtr);
            //printf("%e,%e,%e, %e,%e,%e, %e, %e\n",(ph)->p0, (ph)->p1, (ph)->p2, ph->p3, (ph)->r0, (ph)->r1, (ph)->r2, ph->num_scatt );
            
            (*ph)[i].p0=phHolder->p0;
            (*ph)[i].p1=phHolder->p1;
            (*ph)[i].p2=phHolder->p2;
            (*ph)[i].p3=phHolder->p3;
            (*ph)[i].r0= phHolder->r0; 
            (*ph)[i].r1=phHolder->r1 ;
            (*ph)[i].r2=phHolder->r2; 
            (*ph)[i].num_scatt=phHolder->num_scatt;
            
        }
        
        free(phHolder);
    }
    else
    {
        *framestart+=1; //if the  checkpoint file saved and the program was inturrupted before the frame variable had just increased and before the scatt_frame iteration was saved, add one to the frame start
    }
    
    fclose(fPtr);
}

void readMcPar(char file[200], double *fps, double *theta_jmin, double *theta_j, double *inj_radius, int *frm0,int *last_frm, int *frm2, int *photon_num, double *ph_weight, char *spect, char *restart)
{
    //function to read mc.par file
	FILE *fptr=NULL;
	char buf[100];
	double theta_deg;
	
	//open file
	fptr=fopen(file,"r");
	//read in frames per sec and other variables outlined in main()
	fscanf(fptr, "%lf",fps);
	//printf("%f\n", *fps );
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",frm0);
	//printf("%d\n", *frm0 );
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",last_frm);
	//printf("%d\n", *last_frm );
    
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",frm2);
    *frm2+=*frm0; //frame to go to is what is given in the file plus the starting frame
	//printf("%d\n", *frm2 );
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",photon_num);
	//printf("%d\n", *photon_num );
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%lf",inj_radius);
	//printf("%lf\n", *inj_radius );
	
	fgets(buf, 100,fptr);
	//theta jmin
	fscanf(fptr, "%lf",&theta_deg);
	*theta_jmin=theta_deg*M_PI/180;
	//printf("%f\n", *theta_jmin );
	
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%lf",&theta_deg);
    *theta_j=theta_deg*M_PI/180;
	//printf("%f\n", *theta_j );
	
	fgets(buf, 100,fptr);
    
    fscanf(fptr, "%lf",ph_weight);
    fgets(buf, 100,fptr);
    
    *spect=getc(fptr);
    fgets(buf, 100,fptr);
    //printf("%c\n",*spect);
    
    *restart=getc(fptr);
    
	//close file
	fclose(fptr);
}

void readAndDecimate(char flash_file[200], double r_inj, double **x, double **y, double **szx, double **szy, double **r,\
 double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number)
{
    //function to read in data from FLASH file
    hid_t  file,dset, space;
    herr_t status;
    hsize_t dims[2]={0,0}; //hold dimension size for coordinate data set (mostly interested in dims[0])
    double **vel_x_buffer=NULL, **vel_y_buffer=NULL, **dens_buffer=NULL, **pres_buffer=NULL, **coord_buffer=NULL, **block_sz_buffer=NULL;
    double *velx_unprc=NULL, *vely_unprc=NULL, *dens_unprc=NULL, *pres_unprc=NULL, *x_unprc=NULL, *y_unprc=NULL, *r_unprc=NULL, *szx_unprc=NULL, *szy_unprc=NULL;
    int  i,j,count,x1_count, y1_count, r_count, **node_buffer=NULL, num_nodes=0;
    double x1[8]={-7.0/16,-5.0/16,-3.0/16,-1.0/16,1.0/16,3.0/16,5.0/16,7.0/16};
    
    file = H5Fopen (flash_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    printf(">> mc.py: Reading positional, density, pressure, and velocity information...\n");

    //printf("Reading coord\n");
    dset = H5Dopen (file, "coordinates", H5P_DEFAULT);
    
    //get dimensions of array and save it
    space = H5Dget_space (dset);
    
    H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
    
    /*
     * Allocate array of pointers to rows. OPTIMIZE HERE: INITALIZE ALL THE BUFFERS AT ONCE IN 1 FOR LOOP
     */
    coord_buffer = (double **) malloc (dims[0] * sizeof (double *));
    
    coord_buffer[0] = (double *) malloc (dims[0] * dims[1] * sizeof (double));
    
    block_sz_buffer= (double **) malloc (dims[0] * sizeof (double *));

    block_sz_buffer[0] = (double *) malloc (dims[0] * COORD_DIM1 * sizeof (double));
    
    node_buffer= (int **) malloc (dims[0] * sizeof (int *));
    node_buffer[0] = (int *) malloc (dims[0] * sizeof (int));
    
    vel_x_buffer= (double **) malloc (dims[0] * sizeof (double *));
    vel_x_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
    
    vel_y_buffer= (double **) malloc (dims[0] * sizeof (double *));
    vel_y_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
    
    dens_buffer= (double **) malloc (dims[0] * sizeof (double *));
    dens_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
    
    pres_buffer= (double **) malloc (dims[0] * sizeof (double *));
    pres_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
    
    /*
     * Set the rest of the pointers to rows to the correct addresses.
     */
     for (i=1; i<dims[0]; i++)
     {
         coord_buffer[i] = coord_buffer[0] + i * dims[1];
         block_sz_buffer[i] = block_sz_buffer[0] + i * COORD_DIM1;
         node_buffer[i] = node_buffer[0] + i ;
         vel_x_buffer[i] = vel_x_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
         vel_y_buffer[i] = vel_y_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
         dens_buffer[i] = dens_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
         pres_buffer[i] = pres_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
     }
     
//    for (i=1; i<dims[0]; i++)
//        coord_buffer[i] = coord_buffer[0] + i * dims[1];
    
    //read data such that first column is x and second column is y
    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,coord_buffer[0]);
    
    //close dataset
    status = H5Sclose (space);
    status = H5Dclose (dset);
    
    //printf("Reading block size\n");

    dset = H5Dopen (file, "block size", H5P_DEFAULT);

//    block_sz_buffer= (double **) malloc (dims[0] * sizeof (double *));
//
//    block_sz_buffer[0] = (double *) malloc (dims[0] * COORD_DIM1 * sizeof (double));
//    
//    for (i=1; i<dims[0]; i++){
//        block_sz_buffer[i] = block_sz_buffer[0] + i * COORD_DIM1;
//    }
    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,block_sz_buffer[0]);
    
    // first column of buffer is x and second column is y
    status = H5Dclose (dset);    

    //printf("Reading node type\n");
    dset = H5Dopen (file, "node type", H5P_DEFAULT);

//    node_buffer= (int **) malloc (dims[0] * sizeof (int *));
//    node_buffer[0] = (int *) malloc (dims[0] * sizeof (int));
//
//    for (i=1; i<dims[0]; i++){
//        node_buffer[i] = node_buffer[0] + i ;
//    }
    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,node_buffer[0]);
    status = H5Dclose (dset);

    //printf("Reading velx\n");
    dset = H5Dopen (file, "velx", H5P_DEFAULT);

//   vel_x_buffer= (double **) malloc (dims[0] * sizeof (double *));
//   vel_x_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
//
//   for (i=1; i<dims[0]; i++)
//   {
//        vel_x_buffer[i] = vel_x_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
//   }
   //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,vel_x_buffer[0]);
    status = H5Dclose (dset);

    //printf("Reading vely\n");
    dset = H5Dopen (file, "vely", H5P_DEFAULT);

//     vel_y_buffer= (double **) malloc (dims[0] * sizeof (double *));
//     vel_y_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
//
//     for (i=1; i<dims[0]; i++)
//     {
//        vel_y_buffer[i] = vel_y_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
//     }
    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,vel_y_buffer[0]);
    status = H5Dclose (dset);
    
    //printf("Reading dens\n");
    dset = H5Dopen (file, "dens", H5P_DEFAULT);

//    dens_buffer= (double **) malloc (dims[0] * sizeof (double *));
//    dens_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
//     for (i=1; i<dims[0]; i++)
//     {
//        dens_buffer[i] = dens_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
//     }
     
    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,dens_buffer[0]);
    status = H5Dclose (dset);    
    
    //printf("Reading pres\n");

    dset = H5Dopen (file, "pres", H5P_DEFAULT);

//    pres_buffer= (double **) malloc (dims[0] * sizeof (double *));
//    pres_buffer[0]= (double *) malloc (dims[0] * PROP_DIM1  *PROP_DIM2*PROP_DIM3* sizeof (double));
//     for (i=1; i<dims[0]; i++)
//     {
//        pres_buffer[i] = pres_buffer[0] + i * PROP_DIM1*PROP_DIM2*PROP_DIM3;
//     }
    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,pres_buffer[0]);
    status = H5Dclose (dset);

    status = H5Fclose (file);
    
    printf(">> Selecting good node types (=1)\n");
    //find out how many good nodes there are
    for (i=0;i<dims[0];i++)
    {
        if (node_buffer[i][0]==1 ){
            num_nodes++;
        }
    }
    
    //allocate memory for arrays to hold unprocessed data
    pres_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
        
    dens_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));

    velx_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
     
    vely_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
     
    x_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));

    y_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
    
    r_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
    
    szx_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
    
    szy_unprc=malloc (num_nodes* PROP_DIM1  *PROP_DIM2*PROP_DIM3 * sizeof (double ));
    
    //find where the good values corresponding to the good gones (=1) and save them to the previously allocated pointers which are 1D arrays
    //also create proper x and y arrays and block size arrays
    //and then free up the buffer memory space
    printf(">> Creating and reshaping arrays\n");
    count=0;
    for (i=0;i<dims[0];i++)
    {
        if (node_buffer[i][0]==1 )
        {
            x1_count=0;
            y1_count=0;
            for (j=0;j<(PROP_DIM1*PROP_DIM2*PROP_DIM3);j++)
            {
                *(pres_unprc+count)=pres_buffer[i][j];
                *(dens_unprc+count)=dens_buffer[i][j];
                *(velx_unprc+count)=vel_x_buffer[i][j];
                *(vely_unprc+count)=vel_y_buffer[i][j];
                *(szx_unprc+count)=((block_sz_buffer[i][0])/8)*1e9; //divide by 8 for resolution, multiply by 1e9 to scale properly?
                *(szy_unprc+count)=((block_sz_buffer[i][1])/8)*1e9;
                if (j%8==0)
                {
                    x1_count=0;
                }
                if ((j%8==0) &&  (j!=0))
                {
                    y1_count++;
                }
                *(x_unprc+count)=(coord_buffer[i][0]+block_sz_buffer[i][0]*x1[x1_count])*1e9;
                *(y_unprc+count)=(coord_buffer[i][1]+block_sz_buffer[i][1]*x1[y1_count])*1e9;

                //printf("%d,%d,%d,%d\n",count,j,x1_count,y1_count);
                x1_count++;
                count++;
            }
        }
    }
    free (pres_buffer[0]); free (dens_buffer[0]);free (vel_x_buffer[0]);free (vel_y_buffer[0]); free(coord_buffer[0]);free(block_sz_buffer[0]);free(node_buffer[0]);
    free (pres_buffer);free(dens_buffer);free(vel_x_buffer);free(vel_y_buffer);free(coord_buffer);free(block_sz_buffer);free(node_buffer);
    
    //fill in radius array and find in how many places r > injection radius
    r_count=0;
    for (i=0;i<count;i++)
    {
        *(r_unprc+i)=pow((pow(*(x_unprc+i),2)+pow(*(y_unprc+i),2)),0.5);
        if (*(r_unprc+i)> (0.95*r_inj) )
        {
            r_count++;
        }
    }
        /*
    //find in how many places r > injection radius
    r_count=0;
    for (i=0;i<count;i++)
    {
        if (*(r_unprc+i)> (0.95*r_inj) )
        {
            r_count++;
        }
    }
    */
    //allocate memory to hold processed data
    (*pres)=malloc (r_count * sizeof (double ));
    (*velx)=malloc (r_count * sizeof (double ));
    (*vely)=malloc (r_count * sizeof (double ));
    (*dens)=malloc (r_count * sizeof (double ));
    (*x)=malloc (r_count * sizeof (double ));
    (*y)=malloc (r_count * sizeof (double ));
    (*r)=malloc (r_count * sizeof (double ));
    (*theta)=malloc (r_count * sizeof (double ));
    (*gamma)=malloc (r_count * sizeof (double ));
    (*dens_lab)=malloc (r_count * sizeof (double ));
    (*szx)=malloc (r_count * sizeof (double ));
    (*szy)=malloc (r_count * sizeof (double ));
    (*temp)=malloc (r_count * sizeof (double ));
    
    //assign values based on r> 0.95*r_inj
    j=0;
    for (i=0;i<count;i++)
    {
        if (*(r_unprc+i)> (0.95*r_inj) )
        {
            (*pres)[j]=*(pres_unprc+i);
            (*velx)[j]=*(velx_unprc+i);
            (*vely)[j]=*(vely_unprc+i);
            (*dens)[j]=*(dens_unprc+i);
            (*x)[j]=*(x_unprc+i);
            (*y)[j]=*(y_unprc+i);
            (*r)[j]=*(r_unprc+i);
            (*szx)[j]=*(szx_unprc+i);
            (*szy)[j]=*(szy_unprc+i);
            (*theta)[j]=atan2( *(x_unprc+i) , *(y_unprc+i) );//theta in radians in relation to jet axis
            (*gamma)[j]=pow(pow(1.0-(pow(*(velx_unprc+i),2)+pow(*(vely_unprc+i),2)),0.5),-1); //v is in units of c
            (*dens_lab)[j]= (*(dens_unprc+i)) * (pow(pow(1.0-(pow(*(velx_unprc+i),2)+pow(*(vely_unprc+i),2)),0.5),-1));
            (*temp)[j]=pow(3*(*(pres_unprc+i))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
            j++;
        }
    }
    *number=j;
    
    free(pres_unprc); free(velx_unprc);free(vely_unprc);free(dens_unprc);free(x_unprc); free(y_unprc);free(r_unprc);free(szx_unprc);free(szy_unprc);
    
}


void photonInjection( struct photon **ph, int *ph_num, double r_inj, double ph_weight, char spect, int array_length, double fps, double theta_min, double theta_max,\
double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temps, double *vx, double *vy, gsl_rng * rand)
{
    int i=0, block_cnt=0, *ph_dens=NULL, ph_tot=0, j=0,k=0;
    double ph_dens_calc=0.0, fr_dum=0.0, y_dum=0.0, yfr_dum=0.0, fr_max=0, bb_norm=0, position_phi;
    double com_v_phi, com_v_theta, *p_comv=NULL, *boost=NULL; //comoving phi, theta, comoving 4 momentum for a photon, and boost for photon(to go to lab frame)
    double *l_boost=NULL; //pointer to hold array of lorentz boost, to lab frame, values
    float num_dens_coeff;
    
    if (spect=='w') //from MCRAT paper, w for wien spectrum 
    {
        num_dens_coeff=8.44;
        printf("in wien spectrum\n");
    }
    else
    {
        num_dens_coeff=20.29; //this is for black body spectrum
        //printf("in BB spectrum");
    }
    
    //find how many blocks are near the injection radius within the angles defined in mc.par, get temperatures and calculate number of photons to allocate memory for 
    //and then rcord which blocks have to have "x" amount of photons injected there
    
    for(i=0;i<array_length;i++)
    {
        //look at all boxes in width delta r=c/fps and within angles we are interested in NEED TO IMPLEMENT
            if ((*(r+i) > (r_inj - C_LIGHT/fps))  &&   (*(r+i)  < (r_inj + C_LIGHT/fps)  ) && (*(theta+i)< theta_max) && (*(theta+i) > theta_min) ) 
            {
                block_cnt++;
            }
    }
    //printf("Blocks: %d\n", block_cnt);
    
    //allocate memory to record density of photons for each block
    ph_dens=malloc(block_cnt * sizeof(int));
    
    //calculate the photon density for each block and save it to the array
    j=0;
    ph_tot=0;
    for (i=0;i<array_length;i++)
    {
        //printf("%d\n",i);
        //printf("%e, %e, %e, %e, %e, %e\n", *(r+i),(r_inj - C_LIGHT/fps), (r_inj + C_LIGHT/fps), *(theta+i) , theta_max, theta_min);
            if ((*(r+i) > (r_inj - C_LIGHT/fps))  &&   (*(r+i)  < (r_inj + C_LIGHT/fps)  ) && (*(theta+i)< theta_max) && (*(theta+i) > theta_min) ) 
            {
                ph_dens_calc=num_dens_coeff*2.0*M_PI*(*(x+i))*pow(*(temps+i),3.0)*pow(*(szx+i),2.0) /(ph_weight) ; //a*T^3/(weight) dV, dV=2*PI*x*dx^2, 
                 
                 *(ph_dens+j)=gsl_ran_poisson(rand,ph_dens_calc) ; //choose from poission distribution with mean of ph_dens_calc
                 
                //printf("%d, %lf \n",*(ph_dens+j), ph_dens_calc);
                
                 //sum up all the densities to get total number of photons
                 ph_tot+=(*(ph_dens+j));
                 
                 j++;
            }
    }
        
    //printf("%d\n", ph_tot);
    
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
       if ((*(r+i) > (r_inj - C_LIGHT/fps))  &&   (*(r+i)  < (r_inj + C_LIGHT/fps)  ) && (*(theta+i)< theta_max) && (*(theta+i) > theta_min) )
        {
            
            for(j=0;j<(*(ph_dens+k));j++ )
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
                            fr_max=(C_LIGHT*(*(temps+i)))/(0.29); //max frequency of bb
                            bb_norm=pow((fr_max/(*(temps+i))),2.0)/(exp(PL_CONST*fr_max/K_B/(*(temps+i)))-1); //find value of bb at fr_max
                            yfr_dum=(1.0/bb_norm)*pow((fr_dum/(*(temps+i))),2.0)/(exp((PL_CONST*fr_dum)/(K_B*(*(temps+i)) ))-1);//curve is normalized to vaue of bb @ max frequency
                        }
                        //printf("%lf, %lf,%lf,%e \n",(*(temps+i)),fr_dum, y_dum, yfr_dum);
                        
                    }
                   //printf("%lf\n ",fr_dum);
                   position_phi=gsl_rng_uniform(rand)*2*M_PI;
                   com_v_phi=gsl_rng_uniform(rand)*2*M_PI;
                   com_v_theta=acos((gsl_rng_uniform(rand)*2)-1);
                   //printf("%lf, %lf, %lf\n", position_phi, com_v_phi, com_v_theta);
                   
                   //populate 4 momentum comoving array
                   *(p_comv+0)=PL_CONST*fr_dum/C_LIGHT;
                   *(p_comv+1)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*cos(com_v_phi);
                   *(p_comv+2)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*sin(com_v_phi);
                   *(p_comv+3)=(PL_CONST*fr_dum/C_LIGHT)*cos(com_v_theta);
                   
                    //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
                    *(boost+0)=-1*(*(vx+i))*cos(position_phi);
                    *(boost+1)=-1*(*(vx+i))*sin(position_phi);
                    *(boost+2)=-1*(*(vy+i));
                    //printf("%lf, %lf, %lf\n", *(boost+0), *(boost+1), *(boost+2));
                    
                    //boost to lab frame
                    lorentzBoost(boost, p_comv, l_boost, 'p');
                    //printf("Assignemnt: %e, %e, %e, %e\n", *(l_boost+0), *(l_boost+1), *(l_boost+2),*(l_boost+3));
                   
                (*ph)[ph_tot].p0=(*(l_boost+0));
                (*ph)[ph_tot].p1=(*(l_boost+1));
                (*ph)[ph_tot].p2=(*(l_boost+2));
                (*ph)[ph_tot].p3=(*(l_boost+3));
                (*ph)[ph_tot].r0= (*(x+i))*cos(position_phi); //put photons @ center of box that they are supposed to be in with random phi 
                (*ph)[ph_tot].r1=(*(x+i))*sin(position_phi) ;
                (*ph)[ph_tot].r2=(*(y+i)); //y coordinate in flash becomes z coordinate in MCRaT
                (*ph)[ph_tot].num_scatt=0;
                //printf("%d\n",ph_tot);
                ph_tot++;
            }
            k++;
        }
    }
    
    *ph_num=ph_tot; //save number of photons
    free(ph_dens); free(p_comv);free(boost); free(l_boost);
    
}

void lorentzBoost(double *boost, double *p_ph, double *result, char object)
{
    //function to perform lorentz boost
    //if doing boost for an electron last argument is 'e' and there wont be a check for zero norm
    //if doing boost for a photon  last argument is 'p' and there will be a check for zero norm
    double beta=0, gamma=0, *boosted_p=NULL;
    
    gsl_vector_view b=gsl_vector_view_array(boost, 3); //make boost pointer into vector
    gsl_vector_view p=gsl_vector_view_array(p_ph, 4); //make boost pointer into vector
    gsl_matrix *lambda1= gsl_matrix_calloc (4, 4); //create matrix thats 4x4 to do lorentz boost 
    gsl_vector *p_ph_prime =gsl_vector_calloc(4); //create vestor to hold lorentz boosted vector
    
    //printf("%e, %e, %e, %e\n",gsl_blas_dnrm2(&b.vector), *(boost+0), *(boost+1), *(boost+2));
    //printf("%e, %e, %e, %e\n",*(p_ph+0), *(p_ph+1), *(p_ph+2), *(p_ph+3));
    
    //if magnitude of fluid velocity is != 0 do lorentz boost otherwise dont need to do a boost
    if (gsl_blas_dnrm2(&b.vector) > 0)
    {
        //printf("in If\n");
        beta=sqrt(gsl_blas_dnrm2(&b.vector));
        gamma=1.0/sqrt(1-pow(beta, 2.0));
        
        //initalize matrix values
        gsl_matrix_set(lambda1, 0,0, gamma);
        gsl_matrix_set(lambda1, 0,1,  -1*gsl_vector_get(&b.vector,0)*gamma);
        gsl_matrix_set(lambda1, 0,2,  -1*gsl_vector_get(&b.vector,1)*gamma);
        gsl_matrix_set(lambda1, 0,3,  -1*gsl_vector_get(&b.vector,2)*gamma);
        gsl_matrix_set(lambda1, 1,1,  1+((gamma-1)*pow(gsl_vector_get(&b.vector,0),2.0)/pow(beta,2.0) ) );
        gsl_matrix_set(lambda1, 1,2,  ((gamma-1)*(gsl_vector_get(&b.vector,0)*  gsl_vector_get(&b.vector,1)/pow(beta,2.0) ) ));
        gsl_matrix_set(lambda1, 1,3,  ((gamma-1)*(gsl_vector_get(&b.vector,0)*  gsl_vector_get(&b.vector,2)/pow(beta,2.0) ) ));
        gsl_matrix_set(lambda1, 2,2,  1+((gamma-1)*pow(gsl_vector_get(&b.vector,1),2.0)/pow(beta,2.0) ) );
        gsl_matrix_set(lambda1, 2,3,  ((gamma-1)*(gsl_vector_get(&b.vector,1)*  gsl_vector_get(&b.vector,2)/pow(beta,2.0)) ) );
        gsl_matrix_set(lambda1, 3,3,  1+((gamma-1)*pow(gsl_vector_get(&b.vector,2),2.0)/pow(beta,2.0) ) );
        
        gsl_matrix_set(lambda1, 1,0, gsl_matrix_get(lambda1,0,1));
        gsl_matrix_set(lambda1, 2,0, gsl_matrix_get(lambda1,0,2));
        gsl_matrix_set(lambda1, 3,0, gsl_matrix_get(lambda1,0,3));
        gsl_matrix_set(lambda1, 2,1, gsl_matrix_get(lambda1,1,2));
        gsl_matrix_set(lambda1, 3,1, gsl_matrix_get(lambda1,1,3));
        gsl_matrix_set(lambda1, 3,2, gsl_matrix_get(lambda1,2,3));
        
        gsl_blas_dgemv(CblasNoTrans, 1, lambda1, &p.vector, 0, p_ph_prime );
        //printf("Before Check: %e %e %e %e\n ",gsl_vector_get(p_ph_prime, 0), gsl_vector_get(p_ph_prime, 1), gsl_vector_get(p_ph_prime, 2), gsl_vector_get(p_ph_prime, 3));
        
        //double check vector for 0 norm condition if photon
        if (object == 'p')
        {
            //printf("In if\n");
            boosted_p=zeroNorm(gsl_vector_ptr(p_ph_prime, 0));
        }
        else
        {
            boosted_p=gsl_vector_ptr(p_ph_prime, 0);
        }
        
        //printf("After Check: %e %e %e %e\n ", *(boosted_p+0),*(boosted_p+1),*(boosted_p+2),*(boosted_p+3) );
    }
    else
    {
        //printf("in else");
         //double check vector for 0 norm condition
         if (object=='p')
         {
            boosted_p=zeroNorm(p_ph);
         }
         else
         {
            boosted_p=gsl_vector_ptr(p_ph_prime, 0);
         }
    }
    //assign values to result
    *(result+0)=*(boosted_p+0);
    *(result+1)=*(boosted_p+1);
    *(result+2)=*(boosted_p+2);
    *(result+3)=*(boosted_p+3);
    
    //free up memory
    //free(boosted_p);
    gsl_matrix_free (lambda1); gsl_vector_free(p_ph_prime);
}

double *zeroNorm(double *p_ph)
{
    //ensures zero norm condition of photon 4 monetum is held
    int i=0;
    double normalizing_factor=0;
    gsl_vector_view p=gsl_vector_view_array((p_ph+1), 3); //make last 3 elements of p_ph pointer into vector
    
    if (pow(*(p_ph+0), 2.0) != pow(gsl_blas_dnrm2(&p.vector ),2.0) )
    {
        normalizing_factor=(gsl_blas_dnrm2(&p.vector ));
        //printf("in zero norm if\n");
        //go through and correct 4 momentum assuming the energy is correct
        for (i=1;i<4;i++)
        {
            *(p_ph+i)= ((*(p_ph+i))/(normalizing_factor))*(*(p_ph+0));
        }
        
    }
    /*
     if (pow((*(p_ph+0)),2) != (  pow((*(p_ph+1)),2)+pow((*(p_ph+2)),2)+pow((*(p_ph+3)),2) ) )
        {
            printf("This isnt normalized in the function\nThe difference is: %e\n", pow((*(p_ph+0)),2) - (  pow((*(p_ph+1)),2)+pow((*(p_ph+2)),2)+pow((*(p_ph+3)),2) )  );
        }
    */ //normalized within a factor of 10^-53
    return p_ph;
}

int findNearestPropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *velx,  double *vely, double *dens_lab,\
    double *temp, double *n_dens_lab, double *n_vx, double *n_vy,double *n_temp, gsl_rng * rand)
{
    
    int i=0, j=0, min_index=0;
    double ph_x=0, ph_y=0, ph_phi=0, dist=0, dist_min=1e12;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates
    double ph_v_norm=0, fl_v_norm=0;
    double n_cosangle=0, n_dens_lab_tmp=0,n_vx_tmp=0, n_vy_tmp=0, n_temp_tmp=0 ;
    double rnd_tracker=0, n_dens_lab_min=0, n_vx_min=0, n_vy_min=0, n_temp_min=0;
    
        int index=0;
        double mfp=0,min_mfp=0, beta=0;
       
    //go through each photon and find the blocks around it and then get the distances to all of those blocks and choose the one thats the shortest distance away
    //can optimize here, exchange the for loops and change condition to compare to each of the photons is the radius of the block is .95 (or 1.05) times the min (max) photon radius
    //or just parallelize this part here
    min_mfp=1e12;
    #pragma omp parallel for firstprivate( ph_x, ph_y, ph_phi, dist_min, dist, j, min_index, n_dens_lab_tmp,n_vx_tmp, n_vy_tmp,  n_temp_tmp, fl_v_x, fl_v_y, fl_v_z, fl_v_norm, ph_v_norm, n_cosangle, mfp, beta, rnd_tracker) private(i) shared(min_mfp )
    for (i=0;i<num_ph; i++)
    {
        //printf("%e,%e\n", ((ph+i)->r0), ((ph+i)->r1));
        ph_x=pow(pow(((ph+i)->r0),2.0)+pow(((ph+i)->r1),2.0), 0.5); //convert back to FLASH x coordinate
        ph_y=((ph+i)->r2);
        //printf("ph_x:%e, ph_y:%e\n", ph_x, ph_y);
        ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
        
        dist_min=1e12;//set dist to impossible value to make sure at least first distance calulated is saved
        for(j=0;j<array_num;j++)
        {
            //if the distance between them is within 3e9, to restrict number of possible calculations,  calulate the total distance between the box and photon 
            if ( (fabs(ph_x- (*(x+j)))<3e9) && (fabs(ph_y- (*(y+j)))<3e9))
            {
                //printf("In if statement\n");
                dist= pow(pow(ph_x- (*(x+j)), 2.0) + pow(ph_y- (*(y+j)) , 2.0),0.5);
                //printf("Dist calculated as: %e, index: %d\n", dist, j);
                //printf("In outer if statement, OLD: %e, %d\n", dist_min, min_index);
                
                if((dist<dist_min))
                {
                    //printf("In innermost if statement, OLD: %e, %d\n", dist_min, min_index);
                    dist_min=dist; //save new minimum distance
                    min_index=j; //save index
                    //printf("New Min dist: %e, New min Index: %d\n", dist_min, min_index);
                }
                
            }
        }
        //save values
        /*
        *(n_dens_lab+i)= (*(dens_lab+min_index));
        *(n_vx+i)= (*(velx+min_index));
        *(n_vy+i)= (*(vely+min_index));
        *(n_temp+i)= (*(temp+min_index));
        */
        (n_dens_lab_tmp)= (*(dens_lab+min_index));
        (n_vx_tmp)= (*(velx+min_index));
        (n_vy_tmp)= (*(vely+min_index));
        (n_temp_tmp)= (*(temp+min_index));
        
        
        fl_v_x=(*(velx+min_index))*cos(ph_phi);
        fl_v_y=(*(velx+min_index))*sin(ph_phi);
        fl_v_z=(*(vely+min_index));
        
        fl_v_norm=pow(pow(fl_v_x, 2.0)+pow(fl_v_y, 2.0)+pow(fl_v_z, 2.0), 0.5);
        ph_v_norm=pow(pow(((ph+i)->p1), 2.0)+pow(((ph+i)->p2), 2.0)+pow(((ph+i)->p3), 2.0), 0.5);
        
        //(*(n_cosangle+i))=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //find cosine of the angle between the photon and the fluid velocities via a dot product
        (n_cosangle)=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm );
        
        beta=pow((pow((n_vx_tmp),2)+pow((n_vy_tmp),2)),0.5);
        //put this in to double check that random number is between 0 and 1 (exclusive) because there was a problem with this for parallel case
        rnd_tracker=0;
        while (rnd_tracker<=0 || rnd_tracker>=1)
        {
            rnd_tracker=gsl_rng_uniform_pos(rand);
        }
        mfp=(-1)*(M_P/((n_dens_lab_tmp))/THOMP_X_SECT/(1.0-beta*((n_cosangle))))*log(rnd_tracker) ; //calulate the mfp and then multiply it by the ln of a random number to simulate distribution of mean free paths 
        if (mfp<0)
        {
            printf("\nThread: %d Photon: %d mfp: %e  cos_angle: %e beta: %e dens_lab: %e rnd_tracker: %e\n\n",omp_get_thread_num(), i, mfp, n_cosangle , beta,n_dens_lab_tmp, rnd_tracker );
        }
        
        #pragma omp critical 
        if ( mfp<min_mfp)
        {
            min_mfp=mfp;
            n_dens_lab_min= n_dens_lab_tmp;
            n_vx_min= n_vx_tmp;
            n_vy_min= n_vy_tmp;
            n_temp_min= n_temp_tmp;
            index=i;
            //printf("Thread is %d. new min: %e for photon %d with block properties: %e, %e, %e\n", omp_get_thread_num(), mfp, index, n_vx_tmp, n_vy_tmp, n_temp_tmp);
            #pragma omp flush(min_mfp)
        }

        
    }
    *(n_dens_lab)= n_dens_lab_min;
    *(n_vx)= n_vx_min;
    *(n_vy)= n_vy_min;
    *(n_temp)= n_temp_min;
    (*time_step)=min_mfp/C_LIGHT;
    return index;
    
}

void updatePhotonPosition(struct photon *ph, int num_ph, double t)
{
    //move photons by speed of light
 
    int i=0;
    double old_position=0, new_position=0;
    
    
    
    for (i=0;i<num_ph;i++)
    {
            old_position= pow(  pow(ph->r0,2)+pow(ph->r1,2)+pow(ph->r2,2), 0.5 );
            
            ((ph+i)->r0)+=(((ph+i)->p1)/((ph+i)->p0))*C_LIGHT*t; //update x position
            
            ((ph+i)->r1)+=(((ph+i)->p2)/((ph+i)->p0))*C_LIGHT*t;//update y
            
            ((ph+i)->r2)+=(((ph+i)->p3)/((ph+i)->p0))*C_LIGHT*t;//update z
            
            new_position= pow(  pow(ph->r0,2)+pow(ph->r1,2)+pow(ph->r2,2), 0.5 );
            
            if ((new_position-old_position)/t > C_LIGHT)
            {
                printf("PHOTON NUMBER %d IS SUPERLUMINAL. ITS SPEED IS %e c.\n", i, ((new_position-old_position)/t)/C_LIGHT);
            }
            //printf("In update  function: %e, %e, %e, %e, %e, %e, %e\n",((ph+i)->r0), ((ph+i)->r1), ((ph+i)->r2), t, ((ph+i)->p1)/((ph+i)->p0), ((ph+i)->p2)/((ph+i)->p0), ((ph+i)->p3)/((ph+i)->p0) );  
    }
        
    //printf("In update  function: %e, %e, %e, %e\n",t, ((ph)->p1)/((ph)->p0), ((ph)->p2)/((ph)->p0), ((ph)->p3)/((ph)->p0) );    
    
}


void photonScatter(struct photon *ph, double flash_vx, double flash_vy, double fluid_temp, gsl_rng * rand)
{
    //function to perform single photon scattering
    double ph_phi=0;    
    double *ph_p=malloc(4*sizeof(double)); //pointer to hold only photon 4 momentum @ start
    double *el_p_comov=malloc(4*sizeof(double));//pointer to hold the electron 4 momenta in comoving frame
    double *ph_p_comov=malloc(4*sizeof(double));//pointer to hold the comoving photon 4 momenta
    double *fluid_beta=malloc(3*sizeof(double));//pointer to hold fluid velocity vector
    double *negative_fluid_beta=malloc(3*sizeof(double));//pointer to hold negative fluid velocity vector
    
    /*
    printf("%p\n", ph_p);
    ph_p=calloc(4*sizeof(double),0); 
    printf("%p\n", ph_p);
    el_p_comov=calloc(4*sizeof(double),0);
    printf("%p\n", el_p_comov);
    ph_p_comov=calloc(4*sizeof(double),0);
    printf("%p\n", ph_p_comov);
    fluid_beta=calloc(3*sizeof(double),0);
    printf("%p\n", fluid_beta);
    negative_fluid_beta=calloc(3*sizeof(double),0);
    printf("Done calling calloc\n");
     */
    ph_phi=atan2((ph->r1), ((ph->r0)));
    //printf("ph_phi=%e\n", ph_phi);
    //convert flash coordinated into MCRaT coordinates
    //printf("Getting fluid_beta\n");
    
    (*(fluid_beta+0))=flash_vx*cos(ph_phi);
    (*(fluid_beta+1))=flash_vx*sin(ph_phi);
    (*(fluid_beta+2))=flash_vy;
    
    //fill in photon 4 momentum 
    //printf("filling in 4 momentum in photonScatter\n");
    *(ph_p+0)=(ph->p0);
    *(ph_p+1)=(ph->p1);
    *(ph_p+2)=(ph->p2);
    *(ph_p+3)=(ph->p3);
    
    //first we bring the photon to the fluid's comoving frame
    lorentzBoost(fluid_beta, ph_p, ph_p_comov, 'p');
    //printf("Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
    //printf("New: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
    
    //second we generate a thermal electron at the correct temperature
    singleElectron(el_p_comov, fluid_temp, ph, rand);
    //printf("Outside: %e, %e, %e,%e\n", *(el_p_comov+0), *(el_p_comov+1), *(el_p_comov+2), *(el_p_comov+3));
    
    //third we perform the scattering and save scattered photon 4 monetum in ph_p_comov @ end of function
    singleComptonScatter(el_p_comov, ph_p_comov, rand);
    //printf("Middle: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
    
    //fourth we bring the photon back to the lab frame
     *(negative_fluid_beta+0)=-1*( *(fluid_beta+0));
     *(negative_fluid_beta+1)=-1*( *(fluid_beta+1));
     *(negative_fluid_beta+2)=-1*( *(fluid_beta+2));
    lorentzBoost(negative_fluid_beta, ph_p_comov, ph_p, 'p');
    //printf("Newest: %e, %e, %e,%e\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3));
    //printf("Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
    //printf("Old: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
    
    //assign the photon its new lab 4 momentum
    (ph->p0)=(*(ph_p+0));
    (ph->p1)=(*(ph_p+1));
    (ph->p2)=(*(ph_p+2));
    (ph->p3)=(*(ph_p+3));
    //printf("Done assigning values to original struct\n");

    free(el_p_comov); 
    //printf("done here\n");
    free(ph_p_comov);
    //printf("done here\n");
    free(fluid_beta); // ?maybe not? getting an error checksum for freed object - object was probably modified after being freed.
    //printf("done here\n");
    free(negative_fluid_beta);
    //printf("done here\n");
    free(ph_p);
    //printf("done here\n");
    ph_p=NULL;negative_fluid_beta=NULL;ph_p_comov=NULL; el_p_comov=NULL;
}

void singleElectron(double *el_p, double temp, struct photon *ph, gsl_rng * rand)
{
    //generates an electron with random energy 
    
    double factor=0, gamma=0;
    double y_dum=0, f_x_dum=0, x_dum=0, beta_x_dum=0, beta=0, phi=0, theta=0, ph_theta=0, ph_phi=0;
    gsl_matrix *rot= gsl_matrix_calloc (3, 3); //create matrix thats 3x3 to do rotation 
    gsl_vector_view el_p_prime ; //create vector to hold rotated electron 4 momentum
    gsl_vector *result=gsl_vector_alloc (3);
    
    if (temp>= 1e7)
    {
        //printf("In if\n");
        factor=K_B*temp/(M_EL*pow(C_LIGHT,2.0));
        y_dum=1; //initalize loop to get a random gamma from the distribution of electron velocities
        f_x_dum=0;
        while ((isnan(f_x_dum) !=0) || (y_dum>f_x_dum) )
        {
            
            x_dum=gsl_rng_uniform_pos(rand)*(1+100*factor);
            beta_x_dum=pow(1-(pow(x_dum, -2.0)) ,0.5);
            y_dum=gsl_rng_uniform(rand)/2.0;
            
            f_x_dum=pow(x_dum,2)*(beta_x_dum/gsl_sf_bessel_Kn (2, 1.0/factor))*exp(-1*x_dum/factor); //not sure if this is right is giving small values of gamma -> beta=nan
            //printf("xdum: %e, f_x_dum: %e\n", x_dum, f_x_dum);
        }
        gamma=x_dum;
        
    }
    else
    {

        //printf("In else\n");
        factor=pow(K_B*temp/M_EL,0.5);
        //calculate a random gamma from 3 random velocities drawn from a gaussian distribution with std deviation of "factor"
        gamma=pow( 1- (pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+ pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)  ) ,-0.5);
    }
    
    //printf("%e\n",gamma);
    
    beta=pow( 1- (1/pow(  gamma,2.0 ))  ,0.5);
    //printf("Beta is: %e in singleElectron\n", beta);
    phi=gsl_rng_uniform(rand)*2*M_PI;
    
    y_dum=1; //initalize loop to get a random theta
    f_x_dum=0;
    while (y_dum>f_x_dum)
    {
        y_dum=gsl_rng_uniform(rand)*1.3;
        x_dum=gsl_rng_uniform(rand)*M_PI;
        f_x_dum=sin(x_dum)*(1-(beta*cos(x_dum)));
    }
    theta=x_dum;
    
    //fill in electron 4 momentum NOT SURE WHY THE ORDER IS AS SUCH SEEMS TO BE E/c, pz,py,px!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    *(el_p+0)=gamma*M_EL*C_LIGHT;
    *(el_p+1)=gamma*M_EL*C_LIGHT*beta*cos(theta);
    *(el_p+2)=gamma*M_EL*C_LIGHT*beta*sin(theta)*sin(phi);
    *(el_p+3)=gamma*M_EL*C_LIGHT*beta*sin(theta)*cos(phi);
    
    //printf("Old: %e, %e, %e,%e\n", *(el+0), *(el+1), *(el+2), *(el+3));
    
    el_p_prime=gsl_vector_view_array((el_p+1), 3);
    
    //find angles of photon NOT SURE WHY WERE CHANGING REFERENCE FRAMES HERE???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ph_phi=atan2(ph->p2, ph->p3);
    ph_theta=atan2(pow( pow(ph->p2,2)+  pow(ph->p3,2) , 0.5) ,ph->p1);
    
    //printf("Calculated Photon phi and theta in singleElectron:%lf, %lf\n", ph_phi, ph_theta);
    
    //fill in rotation matrix to rotate around x axis to get rid of phi angle
    gsl_matrix_set(rot, 1,1,1);
    gsl_matrix_set(rot, 2,2,cos(ph_theta));
    gsl_matrix_set(rot, 0,0,cos(ph_theta));
    gsl_matrix_set(rot, 0,2,-sin(ph_theta));
    gsl_matrix_set(rot, 2,0,sin(ph_theta));
    gsl_blas_dgemv(CblasNoTrans, 1, rot, &el_p_prime.vector, 0, result);
    /*
    printf("Middle EL: %e, %e, %e,%e\n", *(el+0), *(el+1), *(el+2), *(el+3));
    printf("Middle: %e, %e, %e,%e\n", *(el+0), gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2));
    printf("Middle EL_P_vec: %e, %e, %e,%e\n", *(el+0), gsl_vector_get(&el_p_prime.vector,0), gsl_vector_get(&el_p_prime.vector,1), gsl_vector_get(&el_p_prime.vector,2));
    */
    gsl_matrix_set_all(rot,0);
    
    gsl_matrix_set(rot, 0,0,1);
    gsl_matrix_set(rot, 1,1,cos(-ph_phi));
    gsl_matrix_set(rot, 2,2,cos(-ph_phi));
    gsl_matrix_set(rot, 1,2,-sin(-ph_phi));
    gsl_matrix_set(rot, 2,1,sin(-ph_phi));
    gsl_blas_dgemv(CblasNoTrans, 1, rot, result, 0, &el_p_prime.vector);
    
    //printf("New: %e, %e, %e,%e\n", *(el+0), *(el+1), *(el+2), *(el+3));
    
    //gsl_rng_free (rand); 
    //printf("freeing pointers in singleElectron\n");
    gsl_matrix_free (rot);gsl_vector_free(result);
    //printf("Done freeing pointers in singleElectron\n");
}

void singleComptonScatter(double *el_comov, double *ph_comov, gsl_rng * rand)
{
    //This routine performs a Compton scattering between a photon and a moving electron.
    int i=0;
    double *el_v=malloc(3*sizeof(double));
    double *negative_el_v=malloc(3*sizeof(double));
    double *ph_p_prime=malloc(4*sizeof(double));//use this to keep track of how the ph 4 momentum changes with each rotation
    double *el_p_prime=malloc(4*sizeof(double));
    double phi0=0, phi1=0, phi=0, theta=0;
    double y_dum, f_x_dum, x_dum;
    gsl_matrix *rot0= gsl_matrix_calloc (3, 3); //create matricies thats 3x3 to do rotations
    gsl_matrix *rot1= gsl_matrix_calloc (3, 3);
    gsl_vector *result0=gsl_vector_alloc (3); //vectors to hold results of rotations
    gsl_vector *result1=gsl_vector_alloc (3); 
    gsl_vector *result=gsl_vector_alloc (4); 
    gsl_vector *whole_ph_p=gsl_vector_alloc (4); 
    gsl_vector_view ph_p ; //create vector to hold comoving photon and electron 4 momentum
    gsl_vector_view el_p ;
    
    //fill in electron velocity array and photon 4 momentum
    *(el_v+0)=(*(el_comov+1))/(*(el_comov+0));
    *(el_v+1)=(*(el_comov+2))/(*(el_comov+0));
    *(el_v+2)=(*(el_comov+3))/(*(el_comov+0));
    
    //lorentz boost into frame where the electron is stationary
    lorentzBoost(el_v, el_comov, el_p_prime, 'e');
    lorentzBoost(el_v, ph_comov, ph_p_prime, 'p');
    
    ph_p=gsl_vector_view_array((ph_p_prime+1), 3);
    el_p=gsl_vector_view_array(el_p_prime,4);
    phi0=atan2(*(ph_p_prime+2), *(ph_p_prime+1) );
    
    //rotate the axes so that the photon incomes along the x-axis
    gsl_matrix_set(rot0, 2,2,1);
    gsl_matrix_set(rot0, 0,0,cos(-phi0));
    gsl_matrix_set(rot0, 1,1,cos(-phi0));
    gsl_matrix_set(rot0, 0,1,-sin(-phi0));
    gsl_matrix_set(rot0, 1,0,sin(-phi0));
    gsl_blas_dgemv(CblasNoTrans, 1, rot0, &ph_p.vector, 0, result0);
    
    //set values of ph_p_prime equal to the result and get new phi from result
    *(ph_p_prime+1)=gsl_vector_get(result0,0);
    *(ph_p_prime+2)=0;//gsl_vector_get(result,1); //just directly setting it to 0 now?
    *(ph_p_prime+3)=gsl_vector_get(result0,2);
    
    phi1=atan2(gsl_vector_get(result0,2), gsl_vector_get(result0,0));
    
    //printf("rotation 1: %e, %e, %e\n",  *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    //printf("rot1: %e, %e, %e,%e\n", *(ph_p_prime+0), gsl_vector_get(&ph_p.vector,0), gsl_vector_get(&ph_p.vector,1), gsl_vector_get(&ph_p.vector,2));
    
    
    //rotate around y to bring it all along x
    gsl_matrix_set(rot1, 1,1,1);
    gsl_matrix_set(rot1, 0,0,cos(-phi1));
    gsl_matrix_set(rot1, 2,2,cos(-phi1));
    gsl_matrix_set(rot1, 0,2,-sin(-phi1));
    gsl_matrix_set(rot1, 2,0,sin(-phi1));
    gsl_blas_dgemv(CblasNoTrans, 1, rot1, &ph_p.vector, 0, result1);
    
    //set values of ph_p_prime equal to the result and get new phi from result
    *(ph_p_prime+1)=*(ph_p_prime+0);//why setting it to the energy?
    *(ph_p_prime+2)=gsl_vector_get(result1,1); 
    *(ph_p_prime+3)=0; //just directly setting it to 0 now?
    
    //printf("rotation 2: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    
    //generate random theta and phi angles for scattering
    phi=gsl_rng_uniform(rand)*2*M_PI;
    
    y_dum=1; //initalize loop to get a random theta
    f_x_dum=0;
    while (y_dum>f_x_dum)
    {
        y_dum=gsl_rng_uniform(rand)*1.09;
        x_dum=gsl_rng_uniform(rand)*M_PI;
        f_x_dum=sin(x_dum)*(1+pow(cos(x_dum),2));
    }
    theta=x_dum;
    
    //perform scattering and compute new 4-momenta of electron and photon
    
    gsl_vector_set(result, 0, (*(ph_p_prime+0))/(1+((*(ph_p_prime+0))/(M_EL*C_LIGHT))*(1-cos(theta)) )); //DOUBLE CHECK HERE!!!!
    gsl_vector_set(result, 1, gsl_vector_get(result,0)*cos(theta) );
    gsl_vector_set(result, 2, gsl_vector_get(result,0)*sin(theta)*sin(phi) );
    gsl_vector_set(result, 3, gsl_vector_get(result,0)*sin(theta)*cos(phi) );
    //printf("%e\n", gsl_vector_get(result,0));
    
    //calculate electron 4 momentum OPTIMIZE HERE: DONT USE A FOR LOOP HERE!!!! Done
    gsl_vector_set(whole_ph_p, 0, (*(ph_p_prime+0)));
    gsl_vector_set(whole_ph_p, 1, (*(ph_p_prime+1)));
    gsl_vector_set(whole_ph_p, 2, (*(ph_p_prime+2)));
    gsl_vector_set(whole_ph_p, 3, (*(ph_p_prime+3)));
    /*
    for (i=0;i<4;i++)
    {
        gsl_vector_set(whole_ph_p, i, (*(ph_p_prime+i)));
    }
    */
     gsl_vector_sub(whole_ph_p,result); //resut is saved into ph_p vector
    gsl_vector_add(&el_p.vector ,whole_ph_p);
    //printf("el_p: %e, %e, %e,%e\n", gsl_vector_get(&el_p.vector,0), gsl_vector_get(&el_p.vector,1), gsl_vector_get(&el_p.vector,2), gsl_vector_get(&el_p.vector,3));
    //printf("ph_p: %e, %e, %e,%e\n", gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2), gsl_vector_get(result,3));
    
    //rotate back to comoving frame
    *(ph_p_prime+0)=gsl_vector_get(result,0);
    *(ph_p_prime+1)=gsl_vector_get(result,1); //set values of photon prime momentum from doing the scattering to use the vector view of it in dot product
    *(ph_p_prime+2)=gsl_vector_get(result,2); 
    *(ph_p_prime+3)=gsl_vector_get(result,3); 
    gsl_matrix_set_all(rot1,0);
    gsl_matrix_set(rot1, 1,1,1);
    gsl_matrix_set(rot1, 0,0,cos(-phi1));
    gsl_matrix_set(rot1, 2,2,cos(-phi1));
    gsl_matrix_set(rot1, 0,2,sin(-phi1));
    gsl_matrix_set(rot1, 2,0,-sin(-phi1));
    gsl_blas_dgemv(CblasNoTrans, 1, rot1, &ph_p.vector, 0, result1);
    
    //set values of ph_p_prime to result1 from undoing 2nd rotation
    *(ph_p_prime+1)=gsl_vector_get(result1,0);
    *(ph_p_prime+2)=gsl_vector_get(result1,1); 
    *(ph_p_prime+3)=gsl_vector_get(result1,2); 
    //printf("Undo rotation 2: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    //ignore the electron, dont care about it, undo the first rotation
    gsl_matrix_set_all(rot0,0);
    gsl_matrix_set(rot0, 2,2,1);
    gsl_matrix_set(rot0, 0,0,cos(-phi0));
    gsl_matrix_set(rot0, 1,1,cos(-phi0));
    gsl_matrix_set(rot0, 0,1,sin(-phi0));
    gsl_matrix_set(rot0, 1,0,-sin(-phi0));
    gsl_blas_dgemv(CblasNoTrans, 1, rot0, &ph_p.vector, 0, result0);
    *(ph_p_prime+1)=gsl_vector_get(result0,0);
    *(ph_p_prime+2)=gsl_vector_get(result0,1); 
    *(ph_p_prime+3)=gsl_vector_get(result0,2); 
    //printf("Undo rotation 1: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    //deboost photon to lab frame
    *(negative_el_v+0)=(-1*(*(el_v+0)));
    *(negative_el_v+1)=(-1*(*(el_v+1)));
    *(negative_el_v+2)=(-1*(*(el_v+2)));
    
    lorentzBoost(negative_el_v, ph_p_prime, ph_comov, 'p');
    //printf("Undo boost 1: %e, %e, %e, %e\n",  *(ph_comov+0), *(ph_comov+1),  *(ph_comov+2),  *(ph_comov+3));
    
    gsl_matrix_free(rot0); gsl_matrix_free(rot1);gsl_vector_free(result0);gsl_vector_free(result1);gsl_vector_free(result);
    //gsl_rng_free (rand); 
    gsl_vector_free(whole_ph_p);free(ph_p_prime);free(el_p_prime);free(el_v); free(negative_el_v);
}


double averagePhotonEnergy(struct photon *ph, int num_ph)
{
    //to calculate average photon energy
    int i=0;
    double sum=0;
    for (i=0;i<num_ph;i++)
    {
        sum+=((ph+i)->p0);
    }
    
    return sum/num_ph;
}

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg )
{
    int temp_max=0, temp_min=-1,  i=0;
    double sum=0;
    
    for (i=0;i<ph_num;i++)
    {
        sum+=((ph+i)->num_scatt);
        
        if (((ph+i)->num_scatt) > temp_max )
        {
            temp_max=((ph+i)->num_scatt);
            //printf("The new max is: %d\n", temp_max);
        }
        
        if ((i==0) || (((ph+i)->num_scatt)<temp_min))
        {
            temp_min=((ph+i)->num_scatt);
            //printf("The new min is: %d\n", temp_min);
        }
    }
    
    *avg=sum/ph_num;
    *max=temp_max;
    *min=temp_min;
    
}
