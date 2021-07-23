//
//  mcrat_flash.c
//  functions that read in Flash hydrodynamic data
//
//  Created by Tyler Parsotan on 7/23/21.
//

#include "mcrat.h"

#define PROP_DIM1 1
#define PROP_DIM2 8
#define PROP_DIM3 8
#define COORD_DIM1 2

void modifyFlashName(char flash_file[STR_BUFFER], char prefix[STR_BUFFER], int frame)
{
    int lim1=0, lim2=0, lim3=0;
    char test[STR_BUFFER]="" ;
    #if DIMENSIONS == 2
    {
        //2D case
        lim1=10;
        lim2=100;
        lim3=1000;
    }
    #else
    {
        //3d case
        lim1=10;
        lim2=100;
        lim3=1000;
    }
    #endif
    
    if (frame<lim1)
    {
        //snprintf(flash_file,sizeof(flash_file), "%s%.3d%d",prefix,000,frame); //FILEPATH,FILEROOT
        snprintf(test,sizeof(test), "%s%s%.3d%d",FILEPATH,FILEROOT,000,frame);
    }
    else if (frame<lim2)
    {
        //snprintf(flash_file,sizeof(flash_file), "%s%.2d%d",prefix,00,frame);
        snprintf(test,sizeof(test), "%s%s%.2d%d",FILEPATH,FILEROOT,00,frame);
    }
    else if (frame<lim3)
    {
        //snprintf(flash_file,sizeof(flash_file), "%s%d%d",prefix,0,frame);
        snprintf(test,sizeof(test), "%s%s%d%d",FILEPATH,FILEROOT,0,frame);
        
    }
    else
    {
        //snprintf(flash_file,sizeof(flash_file), "%s%d",prefix,frame);
        snprintf(test,sizeof(test), "%s%s%d",FILEPATH,FILEROOT,frame);
    }
    strncpy(flash_file, test, sizeof(test));//had to do this workaround for some weird reason
    //printf("test: %s\n", flash_file);
}

void readAndDecimate(char flash_file[STR_BUFFER], struct hydro_dataframe *hydro_data, double r_inj, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr)
{
    //function to read in data from FLASH file
    hid_t  file,dset, space;
    herr_t status;
    hsize_t dims[2]={0,0}; //hold dimension size for coordinate data set (mostly interested in dims[0])
    double **vel_x_buffer=NULL, **vel_y_buffer=NULL, **dens_buffer=NULL, **pres_buffer=NULL, **coord_buffer=NULL, **block_sz_buffer=NULL;
    double *velx_unprc=NULL, *vely_unprc=NULL, *dens_unprc=NULL, *pres_unprc=NULL, *x_unprc=NULL, *y_unprc=NULL, *r_unprc=NULL, *szx_unprc=NULL, *szy_unprc=NULL;
    int  i,j,count,x1_count, y1_count, r_count, **node_buffer=NULL, num_nodes=0, elem_factor=0;
    double x1[8]={-7.0/16,-5.0/16,-3.0/16,-1.0/16,1.0/16,3.0/16,5.0/16,7.0/16};
    double ph_rmin=0, ph_rmax=0, ph_thetamin=0, ph_thetamax=0, r_grid_innercorner=0, r_grid_outercorner=0, theta_grid_innercorner=0, theta_grid_outercorner=0, track_min_r=DBL_MAX, track_max_r=0;
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    

    if (ph_inj_switch==0)
    {
        ph_rmin=min_r;
        ph_rmax=max_r;
        ph_thetamin=min_theta-2*0.017453292519943295; //min_theta - 2*Pi/180 (2 degrees)
        ph_thetamax=max_theta+2*0.017453292519943295; //max_theta + 2*Pi/180 (2 degrees)

    }
    
    file = H5Fopen (flash_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    //ret=H5Pclose(acc_tpl1);
    
    fprintf(fPtr, ">> MCRaT: Reading positional, density, pressure, and velocity information...\n");
    fflush(fPtr);
    //printf("Reading coord\n");
    dset = H5Dopen (file, "coordinates", H5P_DEFAULT);
    
    //get dimensions of array and save it
    space = H5Dget_space (dset);
    
    H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
    
    //status = H5Sclose (space);
    //status = H5Dclose (dset);
    //status = H5Fclose (file);
    
    /*
     * Allocate array of pointers to rows.
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
     

    //read data such that first column is x and second column is y
    //fprintf(fPtr, "Reading Dataset\n");
    //fflush(fPtr);
    //dset = H5Dopen (file, "coordinates", H5P_DEFAULT);
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,coord_buffer[0]);
    
    //close dataset
    status = H5Sclose (space);
    status = H5Dclose (dset);
    
    //printf("Reading block size\n");
    dset = H5Dopen (file, "block size", H5P_DEFAULT);


    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,block_sz_buffer[0]);
    
    // first column of buffer is x and second column is y
    status = H5Dclose (dset);    //status = H5Fclose (file);

    dset = H5Dopen (file, "node type", H5P_DEFAULT);

    status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,node_buffer[0]);
    status = H5Dclose (dset);
    
    
    dset = H5Dopen (file, "velx", H5P_DEFAULT);

   //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,vel_x_buffer[0]);
    status = H5Dclose (dset); //status = H5Fclose (file);

    //printf("Reading vely\n");
    dset = H5Dopen (file, "vely", H5P_DEFAULT);


    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,vel_y_buffer[0]);
    status = H5Dclose (dset); //status = H5Fclose (file);
    
    //printf("Reading dens\n");
    dset = H5Dopen (file, "dens", H5P_DEFAULT);


    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,dens_buffer[0]);
    status = H5Dclose (dset);
    
    //printf("Reading pres\n");
    dset = H5Dopen (file, "pres", H5P_DEFAULT);


    //printf("Reading Dataset\n");
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,pres_buffer[0]);
    status = H5Dclose (dset);
    //H5Pclose(xfer_plist);
    status = H5Fclose (file);

    
    fprintf(fPtr,">> Selecting good node types (=1)\n");
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
    fprintf(fPtr,">> Creating and reshaping arrays\n");
    count=0;
    

    for (i=0;i<dims[0];i++)
    {
        if (node_buffer[i][0]==1 )
        {
            x1_count=0;
            y1_count=0;
            for (j=0;j<(PROP_DIM1*PROP_DIM2*PROP_DIM3);j++)
            {
                *(pres_unprc+count)=pres_buffer[i][j]*HYDRO_P_SCALE;
                *(dens_unprc+count)=dens_buffer[i][j]*HYDRO_D_SCALE;
                *(velx_unprc+count)=vel_x_buffer[i][j];
                *(vely_unprc+count)=vel_y_buffer[i][j];
                *(szx_unprc+count)=((block_sz_buffer[i][0])/8)*HYDRO_L_SCALE; //divide by 8 for resolution, multiply by 1e9 to scale properly?
                *(szy_unprc+count)=((block_sz_buffer[i][1])/8)*HYDRO_L_SCALE;
                if (j%8==0)
                {
                    x1_count=0;
                }
                if ((j%8==0) &&  (j!=0))
                {
                    y1_count++;
                }
                *(x_unprc+count)=(coord_buffer[i][0]+block_sz_buffer[i][0]*x1[x1_count])*HYDRO_L_SCALE;
                *(y_unprc+count)=(coord_buffer[i][1]+block_sz_buffer[i][1]*x1[y1_count])*HYDRO_L_SCALE;

                //printf("%d,%d,%d,%d\n",count,j,x1_count,y1_count);
                x1_count++;
                count++;
            }
        }
    }

    free (pres_buffer[0]); free (dens_buffer[0]);free (vel_x_buffer[0]);free (vel_y_buffer[0]); free(coord_buffer[0]);free(block_sz_buffer[0]);free(node_buffer[0]);
    free (pres_buffer);free(dens_buffer);free(vel_x_buffer);free(vel_y_buffer);free(coord_buffer);free(block_sz_buffer);free(node_buffer);
 

    //fill in radius array and find in how many places r > injection radius
//have single thread execute this while loop and then have inner loop be parallel
    #if CYCLOSYNCHROTRON_SWITCH == ON
        elem_factor=2;
    #else
        elem_factor=0;
    #endif
    r_count=0;
    while (r_count==0)
    {
        r_count=0;
        elem_factor++;
        for (i=0;i<count;i++)
        {
            *(r_unprc+i)=pow((*(x_unprc+i))*(*(x_unprc+i))+(*(y_unprc+i))*(*(y_unprc+i)),0.5);
            
            if (ph_inj_switch==0)
            {
                r_grid_innercorner = pow((*(x_unprc+i) - *(szx_unprc+i)/2.0) * ((*(x_unprc+i) - *(szx_unprc+i)/2.0))+(*(y_unprc+i) - *(szx_unprc+i)/2.0) * (*(y_unprc+i) - *(szx_unprc+i)/2.0),0.5);
                r_grid_outercorner = pow((*(x_unprc+i) + *(szx_unprc+i)/2.0) * ((*(x_unprc+i) + *(szx_unprc+i)/2.0))+(*(y_unprc+i) + *(szx_unprc+i)/2.0) * (*(y_unprc+i) + *(szx_unprc+i)/2.0),0.5);
                
                theta_grid_innercorner = acos( (*(y_unprc+i) - *(szx_unprc+i)/2.0) /r_grid_innercorner); //arccos of y/r for the bottom left corner
                theta_grid_outercorner = acos( (*(y_unprc+i) + *(szx_unprc+i)/2.0) /r_grid_outercorner);
                
                if (((ph_rmin - elem_factor*C_LIGHT/hydro_data->fps) <= r_grid_outercorner) && (r_grid_innercorner  <= (ph_rmax + elem_factor*C_LIGHT/hydro_data->fps) ) && (theta_grid_outercorner >= ph_thetamin) && (theta_grid_innercorner <= ph_thetamax) )
                {
                    r_count++;
                }
            
            }
            else
            {
                if (*(r_unprc+i)> (0.95*r_inj) )
                {
                    r_count++;
                }
            }
        }
        //fprintf(fPtr, "r_count: %d count: %d\n", r_count, count);
    }
    fprintf(fPtr, "Elem factor: %d Ph_rmin: %e rmax: %e Chosen FLASH min_r: %e max_r: %e min_theta: %e degrees max_theta: %e degrees\n", elem_factor, ph_rmin, ph_rmax, ph_rmin - (elem_factor*C_LIGHT/hydro_data->fps), ph_rmax + (elem_factor*C_LIGHT/hydro_data->fps), ph_thetamin*180/M_PI, ph_thetamax*180/M_PI);
    fflush(fPtr);
    
    //allocate memory to hold processed data in the hydro data frame
    (hydro_data->pres)=malloc (r_count * sizeof (double ));
    (hydro_data->dens)=malloc (r_count * sizeof (double ));
    (hydro_data->gamma)=malloc (r_count * sizeof (double ));
    (hydro_data->dens_lab)=malloc (r_count * sizeof (double ));
    (hydro_data->temp)=malloc (r_count * sizeof (double ));
    
    (hydro_data->v0)=malloc (r_count * sizeof (double ));//velx
    (hydro_data->v1)=malloc (r_count * sizeof (double ));//vely
    (hydro_data->r0_size)=malloc (r_count * sizeof (double ));//szx
    (hydro_data->r1_size)=malloc (r_count * sizeof (double ));//szy
    (hydro_data->r0)=malloc (r_count * sizeof (double ));//x
    (hydro_data->r1)=malloc (r_count * sizeof (double ));//y
    (hydro_data->r)=malloc (r_count * sizeof (double ));//r
    (hydro_data->theta)=malloc (r_count * sizeof (double ));//theta


    //assign values based on r> 0.95*r_inj
    j=0;
    for (i=0;i<count;i++)
    {
        if (ph_inj_switch==0)
        {
            r_grid_innercorner = pow((*(x_unprc+i) - *(szx_unprc+i)/2.0) * ((*(x_unprc+i) - *(szx_unprc+i)/2.0))+(*(y_unprc+i) - *(szx_unprc+i)/2.0) * (*(y_unprc+i) - *(szx_unprc+i)/2.0),0.5);
            r_grid_outercorner = pow((*(x_unprc+i) + *(szx_unprc+i)/2.0) * ((*(x_unprc+i) + *(szx_unprc+i)/2.0))+(*(y_unprc+i) + *(szx_unprc+i)/2.0) * (*(y_unprc+i) + *(szx_unprc+i)/2.0),0.5);
            
            theta_grid_innercorner = acos( (*(y_unprc+i) - *(szx_unprc+i)/2.0) /r_grid_innercorner); //arccos of y/r for the bottom left corner
            theta_grid_outercorner = acos( (*(y_unprc+i) + *(szx_unprc+i)/2.0) /r_grid_outercorner);

            if (((ph_rmin - elem_factor*C_LIGHT/hydro_data->fps) <= r_grid_outercorner) && (r_grid_innercorner  <= (ph_rmax + elem_factor*C_LIGHT/hydro_data->fps) ) && (theta_grid_outercorner >= ph_thetamin) && (theta_grid_innercorner <= ph_thetamax))
            {
                ((hydro_data->pres))[j]=*(pres_unprc+i);
                ((hydro_data->v0))[j]=*(velx_unprc+i);
                ((hydro_data->v1))[j]=*(vely_unprc+i);
                
                ((hydro_data->dens))[j]=*(dens_unprc+i);
                ((hydro_data->r0))[j]=*(x_unprc+i);
                ((hydro_data->r1))[j]=*(y_unprc+i);
                ((hydro_data->r))[j]=*(r_unprc+i);
                ((hydro_data->r0_size))[j]=*(szx_unprc+i);
                ((hydro_data->r1_size))[j]=*(szy_unprc+i);
                ((hydro_data->theta))[j]=atan2( *(x_unprc+i) , *(y_unprc+i) );//theta in radians in relation to jet axis
                ((hydro_data->gamma))[j]=pow(pow(1.0-(pow(*(velx_unprc+i),2)+pow(*(vely_unprc+i),2)),0.5),-1); //v is in units of c
                ((hydro_data->dens_lab))[j]= (*(dens_unprc+i)) * (pow(pow(1.0-(pow(*(velx_unprc+i),2)+pow(*(vely_unprc+i),2)),0.5),-1));
                ((hydro_data->temp))[j]=pow(3*(*(pres_unprc+i))/(A_RAD) ,1.0/4.0);
                j++;
                /*
                if (*(r_unprc+i)<track_min_r)
                {
                    track_min_r=*(r_unprc+i);
                }
                
                if (*(r_unprc+i)>track_max_r)
                {
                    track_max_r=*(r_unprc+i);
                }
                 */
            }
        }
        else
        {
            if (*(r_unprc+i)> (0.95*r_inj) )
            {
                ((hydro_data->pres))[j]=*(pres_unprc+i);
                ((hydro_data->v0))[j]=*(velx_unprc+i);
                ((hydro_data->v1))[j]=*(vely_unprc+i);
                ((hydro_data->dens))[j]=*(dens_unprc+i);
                ((hydro_data->r0))[j]=*(x_unprc+i);
                ((hydro_data->r1))[j]=*(y_unprc+i);
                ((hydro_data->r))[j]=*(r_unprc+i);
                ((hydro_data->r0_size))[j]=*(szx_unprc+i);
                ((hydro_data->r1_size))[j]=*(szy_unprc+i);
                ((hydro_data->theta))[j]=atan2( *(x_unprc+i) , *(y_unprc+i) );//theta in radians in relation to jet axis
                ((hydro_data->gamma))[j]=pow(pow(1.0-(pow(*(velx_unprc+i),2)+pow(*(vely_unprc+i),2)),0.5),-1); //v is in units of c
                ((hydro_data->dens_lab))[j]= (*(dens_unprc+i)) * (pow(pow(1.0-(pow(*(velx_unprc+i),2)+pow(*(vely_unprc+i),2)),0.5),-1));
                ((hydro_data->temp))[j]=pow(3*(*(pres_unprc+i))/(A_RAD) ,1.0/4.0);
                j++;
            }
        }
    }
    
    //fprintf(fPtr, "Actual Min and Max Flash grid radii are: %e %e\n", track_min_r, track_max_r);
    //fflush(fPtr);

    hydro_data->num_elements=r_count;

    free(pres_unprc); free(velx_unprc);free(vely_unprc);free(dens_unprc);free(x_unprc); free(y_unprc);free(r_unprc);free(szx_unprc);free(szy_unprc);
    //exit(0);
}
