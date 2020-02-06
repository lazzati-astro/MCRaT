#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <glob.h>
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
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
//#include "mclib_3d.h"
#include "mclib.h"
#include <omp.h>
#include "mpi.h"
#include "mc_synch.h"

#define PROP_DIM1 1
#define PROP_DIM2 8
#define PROP_DIM3 8
#define COORD_DIM1 2
#define R_DIM_2D 9120
#define THETA_DIM_2D 2000

//define constants
const double A_RAD=7.56e-15, C_LIGHT=2.99792458e10, PL_CONST=6.6260755e-27, FINE_STRUCT=7.29735308e-3, CHARGE_EL= 4.8032068e-10;
const double K_B=1.380658e-16, M_P=1.6726231e-24, THOM_X_SECT=6.65246e-25, M_EL=9.1093879e-28 , R_EL=2.817941499892705e-13;

int getOrigNumProcesses(int *counted_cont_procs,  int **proc_array, char dir[200], int angle_rank,  int angle_procs, int last_frame)
{
    int i=0, j=0, val=0, original_num_procs=-1, rand_num=0;
    int frame2=0, framestart=0, scatt_framestart=0, ph_num=0;
    double time=0;
    char mc_chkpt_files[200]="", restrt=""; //define new variable that wont write over the restrt variable in the main part of the code, when its put into the readCheckpoint function
    struct photon *phPtr=NULL; //pointer to array of photons 
    //DIR * dirp;
    //struct dirent * entry;
    //struct stat st = {0};
    glob_t  files;
        
    //if (angle_rank==0)
    {
        //find number of mc_checkpt files there are
        //loop through them and find out which prior processes didnt finish and keep track of which ones didnt
        snprintf(mc_chkpt_files, sizeof(mc_chkpt_files), "%s%s", dir,"mc_chkpt_*" );
        val=glob(mc_chkpt_files, 0, NULL,&files );
    
        //printf("TEST: %s\n", mc_chkpt_files);
    
        //look @ a file by choosing rand int between 0 and files.gl_pathc and if the file exists open and read it to get the actual value for the old number of angle_procs
        srand(angle_rank);
        //printf("NUM_FILES: %d\n",files.gl_pathc);
        
        rand_num=rand() % files.gl_pathc;
        snprintf(mc_chkpt_files, sizeof(mc_chkpt_files), "%s%s%d%s", dir,"mc_chkpt_",  rand_num,".dat" );
        //printf("TEST: %s\n", mc_chkpt_files);
    
        if ( access( mc_chkpt_files, F_OK ) == -1 )
        {
            while(( access( mc_chkpt_files, F_OK ) == -1 ) )
            {
                rand_num=rand() % files.gl_pathc;
                snprintf(mc_chkpt_files, sizeof(mc_chkpt_files), "%s%s%d%s", dir,"mc_chkpt_",  rand_num,".dat" );
                //printf("TEST: %s\n", mc_chkpt_files);
            }
        }
        readCheckpoint(dir, &phPtr, &frame2, &framestart, &scatt_framestart, &ph_num, &restrt, &time, rand_num, &original_num_procs);
    
        //original_num_procs= 70;
        
        
    }
    
    int count_procs[original_num_procs], count=0;
    int cont_procs[original_num_procs];
    //create array of files including any checkpoint file which may not have been created yet b/c old process was still in 1st frame of scattering
    
    for (j=0;j<original_num_procs;j++)
    {
        count_procs[j]=j;
        cont_procs[j]=-1; //set to impossible value for previous mpi process rank that needs to be con't
    }
    
    int limit= (angle_rank != angle_procs-1) ? (angle_rank+1)*original_num_procs/angle_procs : original_num_procs;
    //char mc_chkpt_files[200]="";
    
    printf("Angle ID: %d, start_num: %d, limit: %d\n", angle_rank, (angle_rank*original_num_procs/angle_procs),  limit);
    
    count=0;
    for (j=floor(angle_rank*original_num_procs/angle_procs);j<limit;j++)
    {
        snprintf(mc_chkpt_files, sizeof(mc_chkpt_files), "%s%s%d%s", dir,"mc_chkpt_",  j,".dat" );
        //printf("TEST: %s\n", mc_chkpt_files);
        if ( access( mc_chkpt_files, F_OK ) != -1 )
        {
            readCheckpoint(dir, &phPtr, &frame2, &framestart, &scatt_framestart, &ph_num, &restrt, &time, count_procs[j], &i);
            free(phPtr);
            phPtr=NULL;
            
            if ((framestart<=frame2) && (scatt_framestart<=last_frame)) //add another condition here
            {
                cont_procs[count]=j;
                //printf("ACCEPTED: %s\n", mc_chkpt_files);
                count++;
            }
        }
        else
        {
            cont_procs[count]=j;
            //printf("ACCEPTED: %s\n", mc_chkpt_files);
            count++;
        }
                       
    }
    
    (*proc_array)=malloc (count * sizeof (int )); //allocate space to pointer to hold the old process angle_id's
    count=0;
    for (i=0;i<original_num_procs;i++)
    {
        if (cont_procs[i]!=-1)
        {
            (*proc_array)[count]=cont_procs[i];
            count++;
        }
    }
    
    //save number of old processes this process counted need to be restarted  
    *counted_cont_procs=count;
    
    globfree(& files);
    return original_num_procs;
}


void printPhotons(struct photon *ph, int num_ph, int num_ph_abs, int num_ph_emit, int num_null_ph, int scatt_synch_num_ph, int frame,int frame_inj, int frame_last, char dir[200], int angle_rank, FILE *fPtr )
{
    //function to save the photons' positions and 4 momentum
    
     //now using hdf5 file for each process w/ group structure /(weights or Hydro File #)/(p0,p1,p2,p3, r0, r1, r2, s0, s1, s2, or num_scatt)
    
    //if (STOKES_SWITCH==0), not going to save the polarization info
     
     //open the file if it exists and see if the group exists for the given frame, if frame doesnt exist then write datasets for all photons as extendable
     //if the frame does exist then read information from the prewritten data and then add new data to it as extended chunk
     
     
    int i=0, count=0, rank=1, net_num_ph=num_ph-num_ph_abs-num_null_ph, weight_net_num_ph= num_ph-num_ph_abs-num_null_ph, global_weight_net_num_ph=(frame==frame_inj) ? num_ph-num_ph_abs-num_null_ph : num_ph_emit-num_ph_abs ; //can have more photons absorbed than emitted, weight_net_num_ph=(frame==frame_inj) ? num_ph-num_ph_abs-num_null_ph : scatt_synch_num_ph
    int num_thread=omp_get_num_threads();
    char mc_file[200]="", group[200]="", group_weight[200]="";
    double p0[net_num_ph], p1[net_num_ph], p2[net_num_ph], p3[net_num_ph] , r0[net_num_ph], r1[net_num_ph], r2[net_num_ph], num_scatt[net_num_ph], weight[weight_net_num_ph], global_weight[net_num_ph];
    double s0[net_num_ph], s1[net_num_ph], s2[net_num_ph], s3[net_num_ph], comv_p0[net_num_ph], comv_p1[net_num_ph], comv_p2[net_num_ph], comv_p3[net_num_ph];
    hid_t  file, file_init, dspace, dspace_weight, dspace_global_weight, fspace, mspace, prop, prop_weight, prop_global_weight, group_id;
    hid_t dset_p0, dset_p1, dset_p2, dset_p3, dset_r0, dset_r1, dset_r2, dset_s0, dset_s1, dset_s2, dset_s3, dset_num_scatt, dset_weight, dset_weight_2, dset_comv_p0, dset_comv_p1, dset_comv_p2, dset_comv_p3;
    herr_t status, status_group, status_weight, status_weight_2;
    hsize_t dims[1]={net_num_ph}, dims_weight[1]={weight_net_num_ph}, dims_old[1]={0}; //1 is the number of dimansions for the dataset, called rank

    
    hsize_t maxdims[1]={H5S_UNLIMITED};
    hsize_t      size[1];
    hsize_t      offset[1];
    
    fprintf(fPtr, "num_ph %d num_ph_abs %d num_null_ph %d num_ph_emit %d\nAllocated weight to be %d values large and other arrays to be %d\n",num_ph,num_ph_abs,num_null_ph,num_ph_emit, weight_net_num_ph, net_num_ph);
    
    //save photon data into large arrays, NEED TO KNOW HOW MANY NULL PHOTONS WE HAVE AKA SAVED SPACE THAT AREN'T ACTUALLY PHOTONS TO PROPERLY SAVE SPACE FOR ARRAYS ABOVE
    weight_net_num_ph=0;
    count=0;//used to keep track of weight values since it may not be the same as num_ph
    //#pragma omp parallel for num_threads(num_thread) reduction(+:weight_net_num_ph)
    for (i=0;i<num_ph;i++)
    {
        if ((ph+i)->weight != 0)
        {
            p0[count]= ((ph+i)->p0);
            p1[count]= ((ph+i)->p1);
            p2[count]= ((ph+i)->p2);
            p3[count]= ((ph+i)->p3); 
            r0[count]= ((ph+i)->r0);
            r1[count]= ((ph+i)->r1);
            r2[count]= ((ph+i)->r2);
            s0[count]= ((ph+i)->s0);
            s1[count]= ((ph+i)->s1);
            s2[count]= ((ph+i)->s2);
            s3[count]= ((ph+i)->s3);
            num_scatt[count]= ((ph+i)->num_scatt);
            //if ((frame==frame_inj) || ((scatt_synch_num_ph > 0) && ((ph+i)->type == 'c'))) //if the frame is the same one that the photons were injected in, save the photon weights OR if there are synchrotron photons that havent been absorbed
            {
                weight[weight_net_num_ph]= ((ph+i)->weight);
                weight_net_num_ph++;
                if ((ph+i)->type == 'c')
                {
                    (ph+i)->type = 'o'; //set this to be an old synchrotron scattered photon
                }
            }
            
            if ((frame==frame_last))
            {
                global_weight[count]=((ph+i)->weight);
            }
            
            count++;
        }
    }
    
    
    //make strings for file name and group
    snprintf(mc_file,sizeof(mc_file),"%s%s%d%s",dir,"mc_proc_", angle_rank, ".h5" );
    snprintf(group,sizeof(mc_file),"%d",frame );
    
    //see if file exists, if not create it, if it does just open it
    status = H5Eset_auto(NULL, NULL, NULL); //turn off automatic error printing
    file_init=H5Fcreate(mc_file, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT); //see if the file initially does/doesnt exist
    file=file_init;
    status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr); //turn on auto error printing

    
    if (file_init<0)
    {
        //the file exists, open it with read write 
        file=H5Fopen(mc_file, H5F_ACC_RDWR, H5P_DEFAULT);
        //fprintf(fPtr,"In IF\n");
        
        //see if the group exists
        status = H5Eset_auto(NULL, NULL, NULL);
        status_group = H5Gget_objinfo (file, group, 0, NULL);
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
        
        
        
        /*
        fprintf(fPtr, group);
        if (status_group == 0)
        {   
            fprintf (fPtr, "The group exists.\n");
            //now try to see if there's a weight data set for this group
        }
        else 
        {
            fprintf (fPtr, "The group either does NOT exist\n or some other error occurred.\n"); 
        }
        */
    }
    
    
    if ((file_init>=0) || (status_group != 0) )
    {
        //printf("In IF\n");
        //if the file exists, see if the weight exists
        //snprintf(group_weight,sizeof(group),"/Weight",i );
        status = H5Eset_auto(NULL, NULL, NULL);
        status_weight = H5Gget_objinfo (file, "/Weight", 0, NULL);
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
        
        fprintf(fPtr,"Status of /Weight %d\n", status_weight);
        
        //the file has been newly created or if the group does not exist then  create the group for the frame
        group_id = H5Gcreate2(file, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        
        /* Modify dataset creation properties, i.e. enable chunking  */
        prop = H5Pcreate (H5P_DATASET_CREATE);
        status = H5Pset_chunk (prop, rank, dims);
        
        if ((frame==frame_inj) || (scatt_synch_num_ph > 0))
        {
            prop_weight= H5Pcreate (H5P_DATASET_CREATE);
            status = H5Pset_chunk (prop_weight, rank, dims_weight);
        }
        
        if ((frame==frame_last))
        {
            status = H5Pset_chunk (prop, rank, dims);
        }
    
        /* Create the data space with unlimited dimensions. */
        dspace = H5Screate_simple (rank, dims, maxdims);
        
        dspace_weight=H5Screate_simple (rank, dims_weight, maxdims);
        
        /* Create a new dataset within the file using chunk creation properties.  */
        dset_p0 = H5Dcreate2 (group_id, "P0", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        dset_p1 = H5Dcreate2 (group_id, "P1", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        dset_p2 = H5Dcreate2 (group_id, "P2", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        dset_p3 = H5Dcreate2 (group_id, "P3", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        //if (COMV_SWITCH!=0)
        #if COMV_SWITCH == ON
        {
            dset_comv_p0 = H5Dcreate2 (group_id, "COMV_P0", H5T_NATIVE_DOUBLE, dspace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
            
            dset_comv_p1 = H5Dcreate2 (group_id, "COMV_P1", H5T_NATIVE_DOUBLE, dspace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
            
            dset_comv_p2 = H5Dcreate2 (group_id, "COMV_P2", H5T_NATIVE_DOUBLE, dspace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
            
            dset_comv_p3 = H5Dcreate2 (group_id, "COMV_P3", H5T_NATIVE_DOUBLE, dspace,
                                  H5P_DEFAULT, prop, H5P_DEFAULT);
        }
        #endif
                            
        dset_r0 = H5Dcreate2 (group_id, "R0", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        dset_r1 = H5Dcreate2 (group_id, "R1", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        dset_r2 = H5Dcreate2 (group_id, "R2", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
        
        //if (STOKES_SWITCH!=0)
        #if STOKES_SWITCH == ON
        {
            dset_s0 = H5Dcreate2 (group_id, "S0", H5T_NATIVE_DOUBLE, dspace,
                                H5P_DEFAULT, prop, H5P_DEFAULT);
            
            dset_s1 = H5Dcreate2 (group_id, "S1", H5T_NATIVE_DOUBLE, dspace,
                                H5P_DEFAULT, prop, H5P_DEFAULT);
            
            dset_s2 = H5Dcreate2 (group_id, "S2", H5T_NATIVE_DOUBLE, dspace,
                                H5P_DEFAULT, prop, H5P_DEFAULT);
            
            dset_s3 = H5Dcreate2 (group_id, "S3", H5T_NATIVE_DOUBLE, dspace,
                                H5P_DEFAULT, prop, H5P_DEFAULT);
        }
        #endif
        
        dset_num_scatt = H5Dcreate2 (group_id, "NS", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
                            
        if ((frame==frame_inj) || (scatt_synch_num_ph > 0)) //if the frame is the same one that the photons were injected in, save the photon weights or if there are emitted photons that havent been absorbed
        {
            dset_weight_2 = H5Dcreate2 (group_id, "Weight", H5T_NATIVE_DOUBLE, dspace_weight,
                            H5P_DEFAULT, prop_weight, H5P_DEFAULT); //save the new injected photons' weights
        }
        
        if ((frame==frame_last))
        {
            //if saving the injected photons weight dont have to worry about the major ph_weight thats not in a group
            dset_weight = H5Dcreate2 (file, "Weight", H5T_NATIVE_DOUBLE, dspace,
                                      H5P_DEFAULT, prop, H5P_DEFAULT);
        }
                         
        /* Write data to dataset */
        status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p0);
        
        status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p1);
                        
        status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p2);
                        
        status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p3);
        
        //if (COMV_SWITCH!=0)
        #if COMV_SWITCH == ON
        {
            status = H5Dwrite (dset_comv_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, comv_p0);
            
            status = H5Dwrite (dset_comv_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, comv_p1);
            
            status = H5Dwrite (dset_comv_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, comv_p2);
            
            status = H5Dwrite (dset_comv_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, comv_p3);
        }
        #endif
                        
        status = H5Dwrite (dset_r0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, r0);
        
        status = H5Dwrite (dset_r1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, r1);
                        
        status = H5Dwrite (dset_r2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, r2);
        
        //if (STOKES_SWITCH!=0)
        #if STOKES_SWITCH == ON
        {
            status = H5Dwrite (dset_s0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, s0);
            
            status = H5Dwrite (dset_s1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, s1);
            
            status = H5Dwrite (dset_s2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, s2);
            
            status = H5Dwrite (dset_s3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, s3);
        }
        #endif
        
        status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, num_scatt);
        
        if ((frame==frame_inj) || (scatt_synch_num_ph > 0))
        {
            status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, weight);
            
            status = H5Pclose (prop_weight);
            status = H5Dclose (dset_weight_2);
        }
        
        if ((frame==frame_last))
        {
            printf("Before write\n");
            status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, global_weight);
            printf("After write\n");
        }
        
        status = H5Pclose (prop);
        /*
        if ((status_weight>=0) && (scatt_synch_num_ph > 0) && (frame==frame_last))
        {
            //the /Weight dataset exists (b/c already created it in frame photons were injected in) and we need to do something different to save the emitted synch photons to the dataset
            dset_weight = H5Dopen (file, "Weight", H5P_DEFAULT); //open dataset
            
            //get dimensions of array and save it
            dspace = H5Dget_space (dset_weight);
            
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            
            //extend the dataset
            size[0] = dims_weight[0]+ dims_old[0];
            status = H5Dset_extent (dset_weight, size);
            
           //  Select a hyperslab in extended portion of dataset
            fspace = H5Dget_space (dset_weight);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                          dims_weight, NULL);
            
            // Define memory space
            mspace = H5Screate_simple (rank, dims_weight, NULL);
            
            // Write the data to the extended portion of dataset
            status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, mspace, fspace,
                               H5P_DEFAULT, weight);
            
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
        }
        */
    }
    else
    {
        //if the group already exists then extend it
        //find the size of it now
        
        /* Open an existing group of the specified file. */
        group_id = H5Gopen2(file, group, H5P_DEFAULT);
        dset_p0 = H5Dopen (group_id, "P0", H5P_DEFAULT); //open dataset
    
        //get dimensions of array and save it
        dspace = H5Dget_space (dset_p0);
    
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
        
        //extend the dataset
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_p0, size);
        
        /* Select a hyperslab in extended portion of dataset  */
        fspace = H5Dget_space (dset_p0);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
                                  
        /* Define memory space */
        mspace = H5Screate_simple (rank, dims, NULL);
        
        /* Write the data to the extended portion of dataset  */
        status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, mspace, fspace,
                        H5P_DEFAULT, p0);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
    
        dset_p1 = H5Dopen (group_id, "P1", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_p1);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_p1, size);
        fspace = H5Dget_space (dset_p1);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, p1);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        dset_p2 = H5Dopen (group_id, "P2", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_p2);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_p2, size);
        fspace = H5Dget_space (dset_p2);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, p2);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        dset_p3 = H5Dopen (group_id, "P3", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_p3);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_p3, size);
        fspace = H5Dget_space (dset_p3);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, p3);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        //if (COMV_SWITCH!=0)
        #if COMV_SWITCH == ON
        {
            dset_comv_p0 = H5Dopen (group_id, "COMV_P0", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_comv_p0);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_comv_p0, size);
            fspace = H5Dget_space (dset_comv_p0);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                          dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_comv_p0, H5T_NATIVE_DOUBLE, mspace, fspace,
                               H5P_DEFAULT, comv_p0);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
            dset_comv_p1 = H5Dopen (group_id, "COMV_P1", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_comv_p1);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_comv_p1, size);
            fspace = H5Dget_space (dset_comv_p1);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                          dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_comv_p1, H5T_NATIVE_DOUBLE, mspace, fspace,
                               H5P_DEFAULT, comv_p1);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
            dset_comv_p2 = H5Dopen (group_id, "COMV_P2", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_comv_p2);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_comv_p2, size);
            fspace = H5Dget_space (dset_comv_p2);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                          dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_comv_p2, H5T_NATIVE_DOUBLE, mspace, fspace,
                               H5P_DEFAULT, comv_p2);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
            dset_comv_p3 = H5Dopen (group_id, "COMV_P3", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_comv_p3);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_comv_p3, size);
            fspace = H5Dget_space (dset_comv_p3);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                          dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_comv_p3, H5T_NATIVE_DOUBLE, mspace, fspace,
                               H5P_DEFAULT, comv_p3);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
        }
        #endif
        
        dset_r0 = H5Dopen (group_id, "R0", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_r0);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_r0, size);
        fspace = H5Dget_space (dset_r0);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_r0, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, r0);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        dset_r1 = H5Dopen (group_id, "R1", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_r1);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_r1, size);
        fspace = H5Dget_space (dset_r1);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_r1, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, r1);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        dset_r2 = H5Dopen (group_id, "R2", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_r2);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_r2, size);
        fspace = H5Dget_space (dset_r2);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_r2, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, r2);
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        //if (STOKES_SWITCH!=0)
        #if STOKES_SWITCH == ON
        {
             dset_s0 = H5Dopen (group_id, "S0", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_s0);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_s0, size);
            fspace = H5Dget_space (dset_s0);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_s0, H5T_NATIVE_DOUBLE, mspace, fspace,
                                H5P_DEFAULT, s0);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
            dset_s1 = H5Dopen (group_id, "S1", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_s1);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_s1, size);
            fspace = H5Dget_space (dset_s1);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_s1, H5T_NATIVE_DOUBLE, mspace, fspace,
                                H5P_DEFAULT, s1);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
            dset_s2 = H5Dopen (group_id, "S2", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_s2);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_s2, size);
            fspace = H5Dget_space (dset_s2);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_s2, H5T_NATIVE_DOUBLE, mspace, fspace,
                                H5P_DEFAULT, s2);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
            
            dset_s3 = H5Dopen (group_id, "S3", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_s3);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_s3, size);
            fspace = H5Dget_space (dset_s3);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_s3, H5T_NATIVE_DOUBLE, mspace, fspace,
                                H5P_DEFAULT, s3);
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
        }
        #endif
        
        dset_num_scatt = H5Dopen (group_id, "NS", H5P_DEFAULT); //open dataset
        dspace = H5Dget_space (dset_num_scatt);
        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
        size[0] = dims[0]+ dims_old[0];
        status = H5Dset_extent (dset_num_scatt, size);
        fspace = H5Dget_space (dset_num_scatt);
        offset[0] = dims_old[0];
        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                  dims, NULL); 
        mspace = H5Screate_simple (rank, dims, NULL);
        status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, mspace, fspace,
                            H5P_DEFAULT, num_scatt);
        
        //see if the weights group exists, if it does then we can extend it, otherwise we need to create it and write the new values to it
        snprintf(group_weight,sizeof(group_weight),"Weight",i );
        status = H5Eset_auto(NULL, NULL, NULL);
        status_weight = H5Gget_objinfo (group_id, "Weight", 0, NULL);
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
        
        fprintf(fPtr,"Status of /frame/Weight %d\n", status_weight);
        
        //if (((frame==frame_inj) || (scatt_synch_num_ph > 0)) )
        {
            
            status = H5Sclose (dspace);
            status = H5Sclose (mspace);
            status = H5Sclose (fspace);
    
            if (((frame==frame_last)))
            {
                //make sure to append the newly injected/emitted photons from the most recent set of injected photons to the global weights
            
                dset_weight = H5Dopen (file, "Weight", H5P_DEFAULT); //open dataset
        
                //get dimensions of array and save it
                dspace = H5Dget_space (dset_weight);
        
                status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            
                //extend the dataset
                size[0] = dims_weight[0]+ dims_old[0];
                status = H5Dset_extent (dset_weight, size);
            
                /* Select a hyperslab in extended portion of dataset  */
                fspace = H5Dget_space (dset_weight);
                offset[0] = dims_old[0];
                status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                      dims_weight, NULL);
                
                /* Define memory space */
                mspace = H5Screate_simple (rank, dims_weight, NULL);
            
                /* Write the data to the extended portion of dataset  */
                status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, mspace, fspace,
                                H5P_DEFAULT, weight);
            }
            
            if (status_weight >= 0)
            {
                //will have to create the weight dataset for the new set of phtons that have been injected, although it may already be created since emitting photons now
                //see if the group exists
                status = H5Eset_auto(NULL, NULL, NULL);
                status_weight_2 = H5Gget_objinfo (group_id, "/Weight", 0, NULL);
                status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
                
                if (status_weight_2 < 0)
                {
                    //the dataset doesnt exist
                     /* Modify dataset creation properties, i.e. enable chunking  */
                    prop = H5Pcreate (H5P_DATASET_CREATE);
                    status = H5Pset_chunk (prop, rank, dims);
            
                    /* Create the data space with unlimited dimensions. */
                    dspace = H5Screate_simple (rank, dims, maxdims);
                    
                    dset_weight_2 = H5Dcreate2 (group_id, "Weight", H5T_NATIVE_DOUBLE, dspace,
                                    H5P_DEFAULT, prop, H5P_DEFAULT); //save the new injected photons' weights
                    
                    status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                    H5P_DEFAULT, weight);
                    
                    status = H5Pclose (prop);
                }
                else
                {
                    //it exists and need to modify it
                    dset_weight_2 = H5Dopen (group_id, "Weight", H5P_DEFAULT); //open dataset
                    
                    //get dimensions of array and save it
                    dspace = H5Dget_space (dset_weight_2);
                    
                    status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
                    
                    //extend the dataset
                    size[0] = dims_weight[0]+ dims_old[0];
                    status = H5Dset_extent (dset_weight_2, size);
                    
                    /* Select a hyperslab in extended portion of dataset  */
                    fspace = H5Dget_space (dset_weight_2);
                    offset[0] = dims_old[0];
                    status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                                  dims_weight, NULL);
                    
                    /* Define memory space */
                    mspace = H5Screate_simple (rank, dims_weight, NULL);
                    
                    /* Write the data to the extended portion of dataset  */
                    status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, mspace, fspace,
                                       H5P_DEFAULT, weight);
                }
            }
            else
            {
                fprintf(fPtr, "The frame exists in the hdf5 file but the weight dataset for the frame doesnt exist, therefore creating it.\n");
                fflush(fPtr);
                
                prop_weight= H5Pcreate (H5P_DATASET_CREATE);
                status = H5Pset_chunk (prop_weight, rank, dims_weight);
                
                dspace_weight=H5Screate_simple (rank, dims_weight, maxdims);
                
                dset_weight_2 = H5Dcreate2 (group_id, "Weight", H5T_NATIVE_DOUBLE, dspace_weight,
                                            H5P_DEFAULT, prop_weight, H5P_DEFAULT);
                
                status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                   H5P_DEFAULT, weight);
                
                status = H5Pclose (prop_weight);
            }
            
            status = H5Dclose (dset_weight_2);
        }
        
                        
        
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
    }
    
    
    /* Close resources */
    //status = H5Sclose (dspace);
    status = H5Dclose (dset_p0); status = H5Dclose (dset_p1); status = H5Dclose (dset_p2); status = H5Dclose (dset_p3);
    //if (COMV_SWITCH!=0)
    #if COMV_SWITCH == ON
    {
        status = H5Dclose (dset_comv_p0); status = H5Dclose (dset_comv_p1); status = H5Dclose (dset_comv_p2); status = H5Dclose (dset_comv_p3);
    }
    #endif
    status = H5Dclose (dset_r0); status = H5Dclose (dset_r1); status = H5Dclose (dset_r2);
    //if (STOKES_SWITCH!=0)
    #if STOKES_SWITCH == ON
    {
        status = H5Dclose (dset_s0); status = H5Dclose (dset_s1); status = H5Dclose (dset_s2); status = H5Dclose (dset_s3);
    }
    #endif
    
    status = H5Dclose (dset_num_scatt); 
    //if ((frame==frame_inj) || (scatt_synch_num_ph > 0))
    //{
        //status = H5Dclose (dset_weight_2);
    //} put this right after portion of code that actually uses it
    
    //if ((frame==frame_inj) || ((frame==frame_last) && (status_weight>=0)))
    if ((frame==frame_last))
    {
        status = H5Dclose (dset_weight);
    }
    
    /* Close the group. */
   status = H5Gclose(group_id);
    
    /* Terminate access to the file. */
      status = H5Fclose(file);
    
    //if (status_weight>=0)
    //{
    //    exit(0);
    //}

}

int saveCheckpoint(char dir[200], int frame, int frame2, int scatt_frame, int ph_num,double time_now, struct photon *ph, int last_frame, int angle_rank,int angle_size )
{
    //function to save data necessary to restart simulation if it ends
    //need to save all photon data
    FILE *fPtr=NULL;
    char checkptfile[2000]="";
    char command[2000]="";
    char restart;
    int i=0, success=0;;
    
    
    snprintf(checkptfile,sizeof(checkptfile),"%s%s%d%s",dir,"mc_chkpt_", angle_rank,".dat" );
    
    if ((scatt_frame!=last_frame) && (scatt_frame != frame))
    {
        //quick way to preserve old chkpt file if the new one overwrites the old one and corrupts it for some reason
        snprintf(command, sizeof(command), "%s%s %s_old","exec cp ",checkptfile, checkptfile);
        system(command);
        
        fPtr=fopen(checkptfile, "wb");
        //printf("%s\n", checkptfile);
        
        if (fPtr==NULL)
        {
            printf("Cannot open %s to save checkpoint\n", checkptfile);
            success=1;
        }
        else
        {
            //can call printPhotons here or return an int signifying if the checkpoint save worked
            fwrite(&angle_size, sizeof(int), 1, fPtr);
            restart='c';
            fwrite(&restart, sizeof(char), 1, fPtr);
            //printf("Rank: %d wrote restart %c\n", angle_rank, restart);
            fflush(stdout);
            fwrite(&frame, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote frame\n",  angle_rank);
            fflush(stdout);
            fwrite(&frame2, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote frame2\n",  angle_rank);
            fflush(stdout);
            fwrite(&scatt_frame, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote scatt_frame\n",  angle_rank);
            fflush(stdout);
            fwrite(&time_now, sizeof(double), 1, fPtr);
            //printf("Rank: %d wrote time_now\n",  angle_rank);
            fflush(stdout);
            fwrite(&ph_num, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote ph_num\n",  angle_rank);
            fflush(stdout);
            for(i=0;i<ph_num;i++)
            {
                fwrite((ph+i), sizeof(struct photon ), 1, fPtr);
                //fwrite((ph), sizeof(struct photon )*ph_num, ph_num, fPtr);
            }
            success=0;
        }
        //printf("Rank: %d wrote photons\n",  angle_rank);
        fflush(stdout);
    }
    else if  (scatt_frame == frame)
    {
        snprintf(command, sizeof(command), "%s%s","exec rm ",checkptfile);
        system(command);
        
        fPtr=fopen(checkptfile, "wb");
        //printf("%s\n", checkptfile);
        fflush(stdout);
        
        if (fPtr==NULL)
        {
            printf("Cannot open %s to save checkpoint\n", checkptfile);
            success=1;
        }
        else
        {
            fwrite(&angle_size, sizeof(int), 1, fPtr);
            restart='c';
            fwrite(&restart, sizeof(char), 1, fPtr);
            //printf("Rank: %d wrote restart %c\n", angle_rank, restart);
            fflush(stdout);
            fwrite(&frame, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote frame\n",  angle_rank);
            fflush(stdout);
            fwrite(&frame2, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote frame2\n",  angle_rank);
            fflush(stdout);
            fwrite(&scatt_frame, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote scatt_frame\n",  angle_rank);
            fflush(stdout);
            fwrite(&time_now, sizeof(double), 1, fPtr);
            //printf("Rank: %d wrote time_now\n",  angle_rank);
            fflush(stdout);
            fwrite(&ph_num, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote ph_num\n",  angle_rank);
            fflush(stdout);
            for(i=0;i<ph_num;i++)
            {
                //fwrite((ph), sizeof(struct photon )*ph_num, ph_num, fPtr);
                fwrite((ph+i), sizeof(struct photon ), 1, fPtr);
            }
            //printf("Rank: %d wrote photons\n",  angle_rank);
            success=0;
        }
        fflush(stdout);
        
    }
    else
    {
        //quick way to preserve old chkpt file if the new one overwrites the old one and corrupts it for some reason
        snprintf(command, sizeof(command), "%s%s %s_old","exec cp ",checkptfile, checkptfile);
        system(command);
        
        fPtr=fopen(checkptfile, "wb");
        //printf("%s\n", checkptfile);
        
        if (fPtr==NULL)
        {
            printf("Cannot open %s to save checkpoint\n", checkptfile);
            success=1;
        }
        else
        {
            //just finished last iteration of scatt_frame
            fwrite(&angle_size, sizeof(int), 1, fPtr);
            restart='r';
            fwrite(&restart, sizeof(char), 1, fPtr);
            fwrite(&frame, sizeof(int), 1, fPtr);
            fwrite(&frame2, sizeof(int), 1, fPtr);
            success=0;
        }
    }
    if (success==0)
    {
        fclose(fPtr);
    }
    
    return success;
}

void readCheckpoint(char dir[200], struct photon **ph, int *frame2, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time, int angle_rank, int *angle_size )
{
    //function to read in data from checkpoint file
    FILE *fPtr=NULL;
    char checkptfile[200]="";
    int i=0;
    //int frame, scatt_frame, ph_num, i=0;
    struct photon *phHolder=NULL; //pointer to struct to hold data read in from checkpoint file
    
    
    snprintf(checkptfile,sizeof(checkptfile),"%s%s%d%s",dir,"mc_chkpt_", angle_rank,".dat" );
        
    printf("Checkpoint file: %s\n", checkptfile);
    
    if (access( checkptfile, F_OK ) != -1) //if you can access the file, open and read it
    {
        fPtr=fopen(checkptfile, "rb");
        //if ((angle_rank==2) || (angle_rank==3) || (angle_rank==4) || (angle_rank==5))
        {
            fread(angle_size, sizeof(int), 1, fPtr); //uncomment once I run MCRAT for the sims that didnt save this originally
        }
        fread(restart, sizeof(char), 1, fPtr);
        //printf("%c\n", *restart);
        fread(framestart, sizeof(int), 1, fPtr);
        //printf("%d\n", *framestart);
        fread(frame2, sizeof(int), 1, fPtr);
        
        if((*restart)=='c')
        {
            fread(scatt_framestart, sizeof(int), 1, fPtr);
            
            //if ((riken_switch==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && ((*scatt_framestart)>=3000))
            #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
            if ((*scatt_framestart)>=3000)
            {
                *scatt_framestart+=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1                        
            }
            #else
            {
            *scatt_framestart+=1; //add one to start at the next frame after the simulation was interrrupted
            }
            #endif

            //printf("%d\n", *scatt_framestart);
            fread(time, sizeof(double), 1, fPtr);
            //printf("%e\n", *time);
            fread(ph_num, sizeof(int), 1, fPtr);
            //printf("%d\n", *ph_num);
            
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
                (*ph)[i].comv_p0=phHolder->comv_p0;
                (*ph)[i].comv_p1=phHolder->comv_p1;
                (*ph)[i].comv_p2=phHolder->comv_p2;
                (*ph)[i].comv_p3=phHolder->comv_p3;
                (*ph)[i].r0= phHolder->r0; 
                (*ph)[i].r1=phHolder->r1 ;
                (*ph)[i].r2=phHolder->r2; 
                (*ph)[i].s0=phHolder->s0;  
                (*ph)[i].s1=phHolder->s1;
                (*ph)[i].s2=phHolder->s2;
                (*ph)[i].s3=phHolder->s3;
                (*ph)[i].num_scatt=phHolder->num_scatt;
                (*ph)[i].weight=phHolder->weight;
                (*ph)[i].nearest_block_index= phHolder->nearest_block_index;
                (*ph)[i].type= phHolder->type;
            }
            
            free(phHolder);
        }
        else
        {
            //if ((riken_switch==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && ((*framestart)>=3000))
            #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
            if ((*framestart)>=3000)
            {
                *framestart+=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1                        
            }
            #else
            {
            *framestart+=1; //if the  checkpoint file saved and the program was inturrupted before the frame variable had just increased and before the scatt_frame iteration was saved, add one to the frame start
            }
            #endif
            
            *scatt_framestart=(*framestart);
        }
        
        fclose(fPtr);
    }
    else //if not use default
    {
        //*framestart=(*framestart);
        *scatt_framestart=(*framestart);
        *restart='r';
        
    }
}

void readMcPar(char file[200], double *fluid_domain_x, double *fluid_domain_y, double *fps, double *theta_jmin, double *theta_j, double *d_theta_j, double *inj_radius_small, double *inj_radius_large, int *frm0_small, int *frm0_large, int *last_frm, int *frm2_small,int *frm2_large , double *ph_weight_small,double *ph_weight_large,int *min_photons, int *max_photons, char *spect, char *restart)
{
    //function to read mc.par file
	FILE *fptr=NULL;
	char buf[100]="";
	double theta_deg;
	
	//open file
	fptr=fopen(file,"r");
	//read in frames per sec and other variables outlined in main()
    fscanf(fptr, "%lf",fluid_domain_x);
	//printf("%lf\n", *fluid_domain_x );
	
	fgets(buf, 100,fptr);
    
    fscanf(fptr, "%lf",fluid_domain_y);
	//printf("%lf\n", *fluid_domain_y );
	
	fgets(buf, 100,fptr);
    
	fscanf(fptr, "%lf",fps);
	//printf("%f\n", *fps );
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",frm0_small);
	//printf("%d\n", *frm0_small );
	
	fgets(buf, 100,fptr);
    
    fscanf(fptr, "%d",frm0_large);
	//printf("%d\n", *frm0_large );
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",last_frm);
	//printf("%d\n", *last_frm );
    
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%d",frm2_small);
    *frm2_small+=*frm0_small; //frame to go to is what is given in the file plus the starting frame
	//printf("%d\n", *frm2_small );
	
	fgets(buf, 100,fptr);
	
	//fscanf(fptr, "%d",photon_num); remove photon num because we dont need this
	//printf("%d\n", *photon_num );
    
    fscanf(fptr, "%d",frm2_large);
    *frm2_large+=*frm0_large; //frame to go to is what is given in the file plus the starting frame
    //printf("%d\n", *frm2_large );
	
	fgets(buf, 100,fptr);
	
	//fgets(buf, 100,fptr);
	
	fscanf(fptr, "%lf",inj_radius_small);
	//printf("%lf\n", *inj_radius_small );
	
	fgets(buf, 100,fptr);
    
    fscanf(fptr, "%lf",inj_radius_large);
	//printf("%lf\n", *inj_radius_large );
	
	fgets(buf, 100,fptr);
    
	//theta jmin
	fscanf(fptr, "%lf",&theta_deg);
	*theta_jmin=theta_deg;//*M_PI/180; leave as degrees to manipulate processes 
	//printf("%f\n", *theta_jmin );
	
	
	fgets(buf, 100,fptr);
	
	fscanf(fptr, "%lf",&theta_deg);
    *theta_j=theta_deg;//*M_PI/180;
	//printf("%f\n", *theta_j );
	
	fgets(buf, 100,fptr);
    
    fscanf(fptr, "%lf",d_theta_j);
    //*theta_j=theta_deg;//*M_PI/180;
	//printf("%f\n", *theta_j );
	
	fgets(buf, 100,fptr);
    
    fscanf(fptr, "%lf",ph_weight_small);
    //printf("%f\n", *ph_weight_small );
    fgets(buf, 100,fptr);
    
    fscanf(fptr, "%lf",ph_weight_large);
    fgets(buf, 100,fptr);
    
    fscanf(fptr, "%d",min_photons);
    fgets(buf, 100,fptr);
    
    fscanf(fptr, "%d",max_photons);
    fgets(buf, 100,fptr);
    
    *spect=getc(fptr);
    fgets(buf, 100,fptr);
    //printf("%c\n",*spect);
    
    *restart=getc(fptr);
    fgets(buf, 100,fptr);
    
    //dont need this line fo code for MPI 
    //fscanf(fptr, "%d",num_threads);
    //printf("MAKE SURE THERE IS NO NUM_THREADS LINE IN THE MC.PAR FILE.\n");
    //fgets(buf, 100,fptr);
    
    //fscanf(fptr, "%d",dim_switch);
    //printf("MAKE SURE THERE IS NO DIM_SWITCH LINE IN THE MC.PAR FILE.\n");
    //printf("%d\n",*dim_switch);
    
	//close file
	fclose(fptr);
}

void readAndDecimate(char flash_file[200], double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
 double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr)
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
    int num_thread=omp_get_num_threads();
    
    
    

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
                *(pres_unprc+count)=pres_buffer[i][j];
                *(dens_unprc+count)=dens_buffer[i][j];
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
    elem_factor=0;
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
                
                if (((ph_rmin - elem_factor*C_LIGHT/fps) <= r_grid_outercorner) && (r_grid_innercorner  <= (ph_rmax + elem_factor*C_LIGHT/fps) ) && (theta_grid_outercorner >= ph_thetamin) && (theta_grid_innercorner <= ph_thetamax) )
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
    fprintf(fPtr, "Elem factor: %d Ph_rmin: %e rmax: %e Chosen FLASH min_r: %e max_r: %e min_theta: %e degrees max_theta: %e degrees\n", elem_factor, ph_rmin, ph_rmax, ph_rmin - (elem_factor*C_LIGHT/fps), ph_rmax + (elem_factor*C_LIGHT/fps), ph_thetamin*180/M_PI, ph_thetamax*180/M_PI);
    fflush(fPtr);
    
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
        if (ph_inj_switch==0)
        {
            r_grid_innercorner = pow((*(x_unprc+i) - *(szx_unprc+i)/2.0) * ((*(x_unprc+i) - *(szx_unprc+i)/2.0))+(*(y_unprc+i) - *(szx_unprc+i)/2.0) * (*(y_unprc+i) - *(szx_unprc+i)/2.0),0.5);
            r_grid_outercorner = pow((*(x_unprc+i) + *(szx_unprc+i)/2.0) * ((*(x_unprc+i) + *(szx_unprc+i)/2.0))+(*(y_unprc+i) + *(szx_unprc+i)/2.0) * (*(y_unprc+i) + *(szx_unprc+i)/2.0),0.5);
            
            theta_grid_innercorner = acos( (*(y_unprc+i) - *(szx_unprc+i)/2.0) /r_grid_innercorner); //arccos of y/r for the bottom left corner
            theta_grid_outercorner = acos( (*(y_unprc+i) + *(szx_unprc+i)/2.0) /r_grid_outercorner);

            if (((ph_rmin - elem_factor*C_LIGHT/fps) <= r_grid_outercorner) && (r_grid_innercorner  <= (ph_rmax + elem_factor*C_LIGHT/fps) ) && (theta_grid_outercorner >= ph_thetamin) && (theta_grid_innercorner <= ph_thetamax))
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
    }
    
    //fprintf(fPtr, "Actual Min and Max Flash grid radii are: %e %e\n", track_min_r, track_max_r);
    //fflush(fPtr);

    *number=r_count;

    free(pres_unprc); free(velx_unprc);free(vely_unprc);free(dens_unprc);free(x_unprc); free(y_unprc);free(r_unprc);free(szx_unprc);free(szy_unprc);
    //exit(0);
}


void photonInjection( struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, int array_length, double fps, double theta_min, double theta_max,\
double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temps, double *vx, double *vy, gsl_rng * rand, FILE *fPtr)
{
    int i=0, block_cnt=0, *ph_dens=NULL, ph_tot=0, j=0,k=0;
    double ph_dens_calc=0.0, fr_dum=0.0, y_dum=0.0, yfr_dum=0.0, fr_max=0, bb_norm=0, position_phi, ph_weight_adjusted, rmin, rmax;
    double com_v_phi, com_v_theta, *p_comv=NULL, *boost=NULL; //comoving phi, theta, comoving 4 momentum for a photon, and boost for photon(to go to lab frame)
    double *l_boost=NULL; //pointer to hold array of lorentz boost, to lab frame, values
    float num_dens_coeff;
    double r_grid_innercorner=0, r_grid_outercorner=0, theta_grid_innercorner=0, theta_grid_outercorner=0;
    double position_rand=0, position2_rand=0;
    
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
    
    rmin=r_inj - 0.5*C_LIGHT/fps;
    rmax=r_inj + 0.5*C_LIGHT/fps;
    
    for(i=0;i<array_length;i++)
    {
        #if GEOMETRY == CARTESIAN
            r_grid_innercorner = pow((*(x+i) - *(szx+i)/2.0) * ((*(x+i) - *(szx+i)/2.0))+(*(y+i) - *(szy+i)/2.0) * (*(y+i) - *(szy+i)/2.0),0.5);
            r_grid_outercorner = pow((*(x+i) + *(szx+i)/2.0) * ((*(x+i) + *(szx+i)/2.0))+(*(y+i) + *(szy+i)/2.0) * (*(y+i) + *(szy+i)/2.0),0.5);
            
            theta_grid_innercorner = acos( (*(y+i) - *(szx+i)/2.0) /r_grid_innercorner); //arccos of y/r for the bottom left corner
            theta_grid_outercorner = acos( (*(y+i) + *(szx+i)/2.0) /r_grid_outercorner);
        #elif GEOMETRY == SPHERICAL
            r_grid_innercorner = (*(r+i)) - 0.5 * (*(szx+i));
            r_grid_outercorner = (*(r+i)) + 0.5 * (*(szx+i));
            
            theta_grid_innercorner = (*(theta+i)) - 0.5 * (*(szy+i));
            theta_grid_outercorner = (*(theta+i)) + 0.5 * (*(szy+i));
        #endif

        //look at all boxes in width delta r=c/fps and within angles we are interested in NEED TO IMPLEMENT
        //if ((*(r+i) >= rmin)  &&   (*(r+i)  < rmax  ) && (*(theta+i)< theta_max) && (*(theta+i) >=theta_min) )
        if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max))
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
    ph_weight_adjusted=ph_weight;
    //printf("%d %d\n", max_photons, min_photons);
    while ((ph_tot>max_photons) || (ph_tot<min_photons) )
    {
        j=0;
        ph_tot=0;
        
        for (i=0;i<array_length;i++)
        {
            //printf("%d\n",i);
            //printf("%e, %e, %e, %e, %e, %e\n", *(r+i),(r_inj - C_LIGHT/fps), (r_inj + C_LIGHT/fps), *(theta+i) , theta_max, theta_min);
            #if GEOMETRY == CARTESIAN
                r_grid_innercorner = pow((*(x+i) - *(szx+i)/2.0) * ((*(x+i) - *(szx+i)/2.0))+(*(y+i) - *(szy+i)/2.0) * (*(y+i) - *(szy+i)/2.0),0.5);
                r_grid_outercorner = pow((*(x+i) + *(szx+i)/2.0) * ((*(x+i) + *(szx+i)/2.0))+(*(y+i) + *(szy+i)/2.0) * (*(y+i) + *(szy+i)/2.0),0.5);
                
                theta_grid_innercorner = acos( (*(y+i) - *(szx+i)/2.0) /r_grid_innercorner); //arccos of y/r for the bottom left corner
                theta_grid_outercorner = acos( (*(y+i) + *(szx+i)/2.0) /r_grid_outercorner);
            #elif GEOMETRY == SPHERICAL
                r_grid_innercorner = (*(r+i)) - 0.5 * (*(szx+i));
                r_grid_outercorner = (*(r+i)) + 0.5 * (*(szx+i));
                
                theta_grid_innercorner = (*(theta+i)) - 0.5 * (*(szy+i));
                theta_grid_outercorner = (*(theta+i)) + 0.5 * (*(szy+i));
            #endif

            //if ((*(r+i) >= rmin)  &&   (*(r+i)  < rmax  ) && (*(theta+i)< theta_max) && (*(theta+i) >=theta_min) )
            if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max))
                {
                    #if GEOMETRY == SPHERICAL
                    {
                        ph_dens_calc=(num_dens_coeff*2.0*M_PI*pow(*(r+i),2)*sin(*(theta+i))*pow(*(temps+i),3.0)*(*(szx+i))*(*(szy+i)) /(ph_weight_adjusted))*pow(pow(1.0-(pow(*(vx+i),2)+pow(*(vy+i),2)),0.5),-1); //dV=2 *pi* r^2 Sin(theta) dr dtheta
                    }
                    #else
                    {
                        //using FLASH
                        ph_dens_calc=(4.0/3.0)*(num_dens_coeff*2.0*M_PI*(*(x+i))*pow(*(temps+i),3.0)*(*(szx+i))*(*(szy+i)) /(ph_weight_adjusted))*pow(pow(1.0-(pow(*(vx+i),2)+pow(*(vy+i),2)),0.5),-1) ; //a*T^3/(weight) dV, dV=2*PI*x*dx^2,
                    }
                    #endif
                    
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
            
        }
        else if (ph_tot<min_photons)
        {
            ph_weight_adjusted*=0.5;
            
        }
        
        //printf("dens: %d, photons: %d\n", *(ph_dens+(j-1)), ph_tot);
         
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
        #if GEOMETRY == CARTESIAN
            r_grid_innercorner = pow((*(x+i) - *(szx+i)/2.0) * ((*(x+i) - *(szx+i)/2.0))+(*(y+i) - *(szy+i)/2.0) * (*(y+i) - *(szy+i)/2.0),0.5);
            r_grid_outercorner = pow((*(x+i) + *(szx+i)/2.0) * ((*(x+i) + *(szx+i)/2.0))+(*(y+i) + *(szy+i)/2.0) * (*(y+i) + *(szy+i)/2.0),0.5);
            
            theta_grid_innercorner = acos( (*(y+i) - *(szx+i)/2.0) /r_grid_innercorner); //arccos of y/r for the bottom left corner
            theta_grid_outercorner = acos( (*(y+i) + *(szx+i)/2.0) /r_grid_outercorner);
        #elif GEOMETRY == SPHERICAL
            r_grid_innercorner = (*(r+i)) - 0.5 * (*(szx+i));
            r_grid_outercorner = (*(r+i)) + 0.5 * (*(szx+i));
            
            theta_grid_innercorner = (*(theta+i)) - 0.5 * (*(szy+i));
            theta_grid_outercorner = (*(theta+i)) + 0.5 * (*(szy+i));
        #endif
        
        //if ((*(r+i) >= rmin)  &&   (*(r+i)  < rmax  ) && (*(theta+i)< theta_max) && (*(theta+i) >=theta_min) )
        if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max))
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
                    //printf("i: %d freq:%lf\n ",ph_tot, fr_dum);
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
                    
                    //boost to lab frame
                    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
                    //printf("Assignemnt: %e, %e, %e, %e\n", *(l_boost+0), *(l_boost+1), *(l_boost+2),*(l_boost+3));
                   
                (*ph)[ph_tot].p0=(*(l_boost+0));
                (*ph)[ph_tot].p1=(*(l_boost+1));
                (*ph)[ph_tot].p2=(*(l_boost+2));
                (*ph)[ph_tot].p3=(*(l_boost+3));
                (*ph)[ph_tot].comv_p0=(*(p_comv+0));
                (*ph)[ph_tot].comv_p1=(*(p_comv+1));
                (*ph)[ph_tot].comv_p2=(*(p_comv+2));
                (*ph)[ph_tot].comv_p3=(*(p_comv+3));
                
                //place photons in rand positions within fluid element
                #if GEOMETRY == CARTESIAN
                    position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
                    (*ph)[ph_tot].r0= (*(x+i)+position_rand)*cos(position_phi); //put photons @ center of box that they are supposed to be in with random phi
                    (*ph)[ph_tot].r1=(*(x+i)+position_rand)*sin(position_phi) ;
                    position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0;
                    (*ph)[ph_tot].r2=(*(y+i)+position_rand); //y coordinate in flash becomes z coordinate in MCRaT
                #elif GEOMETRY == SPHERICAL
                    position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
                    position2_rand=gsl_rng_uniform_pos(rand)*(*(szy+i))-(*(szy+i))/2.0;
                
                    (*ph)[ph_tot].r0= (*(r+i)+position_rand)*sin(*(theta+i)+position2_rand)*cos(position_phi); //put photons @ center of box that they are supposed to be in with random phi
                    (*ph)[ph_tot].r1=(*(r+i)+position_rand)*sin(*(theta+i)+position2_rand)*sin(position_phi) ;
                    
                    (*ph)[ph_tot].r2=(*(r+i)+position_rand)*cos(*(theta+i)+position2_rand); //y coordinate in flash becomes z coordinate in MCRaT
                #endif
                
                (*ph)[ph_tot].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1 
                (*ph)[ph_tot].s1=0;
                (*ph)[ph_tot].s2=0;
                (*ph)[ph_tot].s3=0;
                (*ph)[ph_tot].num_scatt=0;
                (*ph)[ph_tot].weight=ph_weight_adjusted;
                (*ph)[ph_tot].nearest_block_index=0;
                (*ph)[ph_tot].type='i'; //i for injected
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

void lorentzBoost(double *boost, double *p_ph, double *result, char object,  FILE *fPtr)
{
    //function to perform lorentz boost
    //if doing boost for an electron last argument is 'e' and there wont be a check for zero norm
    //if doing boost for a photon  last argument is 'p' and there will be a check for zero norm
    double beta=0, gamma=0, *boosted_p=NULL;
    
    gsl_vector_view b=gsl_vector_view_array(boost, 3); //make boost pointer into vector
    gsl_vector_view p=gsl_vector_view_array(p_ph, 4); //make boost pointer into vector
    gsl_matrix *lambda1= gsl_matrix_calloc (4, 4); //create matrix thats 4x4 to do lorentz boost 
    gsl_vector *p_ph_prime =gsl_vector_calloc(4); //create vestor to hold lorentz boosted vector
    
    /*
    fprintf(fPtr,"Boost: %e, %e, %e, %e\n",gsl_blas_dnrm2(&b.vector), *(boost+0), *(boost+1), *(boost+2));
    fflush(fPtr);
    fprintf(fPtr,"4 Momentum to Boost: %e, %e, %e, %e\n",*(p_ph+0), *(p_ph+1), *(p_ph+2), *(p_ph+3));
    fflush(fPtr);
    */
    
    //if magnitude of fluid velocity is != 0 do lorentz boost otherwise dont need to do a boost
    if (gsl_blas_dnrm2(&b.vector) > 0)
    {
        //fprintf(fPtr,"in If\n");
        //fflush(fPtr);
        beta=gsl_blas_dnrm2(&b.vector);
        gamma=1.0/sqrt(1-pow(beta, 2.0));
        //fprintf(fPtr,"Beta: %e\tGamma: %e\n",beta,gamma );
        //fflush(fPtr);
        
        //initalize matrix values
        gsl_matrix_set(lambda1, 0,0, gamma);
        gsl_matrix_set(lambda1, 0,1,  -1*gsl_vector_get(&b.vector,0)*gamma);
        gsl_matrix_set(lambda1, 0,2,  -1*gsl_vector_get(&b.vector,1)*gamma);
        gsl_matrix_set(lambda1, 0,3,  -1*gsl_vector_get(&b.vector,2)*gamma);
        gsl_matrix_set(lambda1, 1,1,  1+((gamma-1)*(gsl_vector_get(&b.vector,0)*gsl_vector_get(&b.vector,0))/(beta*beta) ) );
        gsl_matrix_set(lambda1, 1,2,  ((gamma-1)*(gsl_vector_get(&b.vector,0)*  gsl_vector_get(&b.vector,1)/(beta*beta) ) ));
        gsl_matrix_set(lambda1, 1,3,  ((gamma-1)*(gsl_vector_get(&b.vector,0)*  gsl_vector_get(&b.vector,2)/(beta*beta) ) ));
        gsl_matrix_set(lambda1, 2,2,  1+((gamma-1)*(gsl_vector_get(&b.vector,1)*gsl_vector_get(&b.vector,1))/(beta*beta) ) );
        gsl_matrix_set(lambda1, 2,3,  ((gamma-1)*(gsl_vector_get(&b.vector,1)*  gsl_vector_get(&b.vector,2))/(beta*beta) ) );
        gsl_matrix_set(lambda1, 3,3,  1+((gamma-1)*(gsl_vector_get(&b.vector,2)*gsl_vector_get(&b.vector,2))/(beta*beta) ) );
        
        gsl_matrix_set(lambda1, 1,0, gsl_matrix_get(lambda1,0,1));
        gsl_matrix_set(lambda1, 2,0, gsl_matrix_get(lambda1,0,2));
        gsl_matrix_set(lambda1, 3,0, gsl_matrix_get(lambda1,0,3));
        gsl_matrix_set(lambda1, 2,1, gsl_matrix_get(lambda1,1,2));
        gsl_matrix_set(lambda1, 3,1, gsl_matrix_get(lambda1,1,3));
        gsl_matrix_set(lambda1, 3,2, gsl_matrix_get(lambda1,2,3));
        
        gsl_blas_dgemv(CblasNoTrans, 1, lambda1, &p.vector, 0, p_ph_prime );
        
        /*
        fprintf(fPtr,"Lorentz Boost Matrix 0: %e,%e, %e, %e\n", gsl_matrix_get(lambda1, 0,0), gsl_matrix_get(lambda1, 0,1), gsl_matrix_get(lambda1, 0,2), gsl_matrix_get(lambda1, 0,3));
        fflush(fPtr);
        fprintf(fPtr,"Lorentz Boost Matrix 1: %e,%e, %e, %e\n", gsl_matrix_get(lambda1, 1,0), gsl_matrix_get(lambda1, 1,1), gsl_matrix_get(lambda1, 1,2), gsl_matrix_get(lambda1, 1,3));
        fflush(fPtr);
        fprintf(fPtr,"Lorentz Boost Matrix 2: %e,%e, %e, %e\n", gsl_matrix_get(lambda1, 2,0), gsl_matrix_get(lambda1, 2,1), gsl_matrix_get(lambda1, 2,2), gsl_matrix_get(lambda1, 2,3));
        fflush(fPtr);
        fprintf(fPtr,"Lorentz Boost Matrix 3: %e,%e, %e, %e\n", gsl_matrix_get(lambda1, 3,0), gsl_matrix_get(lambda1, 3,1), gsl_matrix_get(lambda1, 3,2), gsl_matrix_get(lambda1, 3,3));
        fflush(fPtr);
        
        fprintf(fPtr,"Before Check: %e %e %e %e\n ",gsl_vector_get(p_ph_prime, 0), gsl_vector_get(p_ph_prime, 1), gsl_vector_get(p_ph_prime, 2), gsl_vector_get(p_ph_prime, 3));
        fflush(fPtr);
        */
        
        //double check vector for 0 norm condition if photon
        if (object == 'p')
        {
            //fprintf(fPtr,"In if\n");
            boosted_p=zeroNorm(gsl_vector_ptr(p_ph_prime, 0));
        }
        else
        {
            boosted_p=gsl_vector_ptr(p_ph_prime, 0);
        }
        /*
        fprintf(fPtr,"After Check: %e %e %e %e\n ", *(boosted_p+0),*(boosted_p+1),*(boosted_p+2),*(boosted_p+3) );
        fflush(fPtr);
         * */
    }
    else
    {
        /*
        fprintf(fPtr,"in else");
        fflush(fPtr);
         * */
         //double check vector for 0 norm condition
         if (object=='p')
         {
            boosted_p=zeroNorm(p_ph);
         }
         else
         {
             //if 4 momentum isnt for photon and there is no boost to be done, we dont care about normality and just want back what was passed to lorentz boost
            boosted_p=gsl_vector_ptr(&p.vector, 0);
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
    
    if (*(p_ph+0) != gsl_blas_dnrm2(&p.vector ) )
    {
        normalizing_factor=(gsl_blas_dnrm2(&p.vector ));
        //fprintf(fPtr,"in zero norm if\n");
        //fflush(fPtr);
        //go through and correct 4 momentum assuming the energy is correct
        
        *(p_ph+1)= ((*(p_ph+1))/(normalizing_factor))*(*(p_ph+0));
        *(p_ph+2)= ((*(p_ph+2))/(normalizing_factor))*(*(p_ph+0));
        *(p_ph+3)= ((*(p_ph+3))/(normalizing_factor))*(*(p_ph+0));
        
    }
    /*
     if (pow((*(p_ph+0)),2) != (  pow((*(p_ph+1)),2)+pow((*(p_ph+2)),2)+pow((*(p_ph+3)),2) ) )
        {
            printf("This isnt normalized in the function\nThe difference is: %e\n", pow((*(p_ph+0)),2) - (  pow((*(p_ph+1)),2)+pow((*(p_ph+2)),2)+pow((*(p_ph+3)),2) )  );
        }
    */ //normalized within a factor of 10^-53
    return p_ph;
}

int findNearestBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z)
{
    double dist=0, dist_min=1e15, block_dist=0;
    int min_index=0, j=0;
    
    dist_min=1e15;//set dist to impossible value to make sure at least first distance calulated is saved
    block_dist=3e9;
    while (dist_min==1e15) //if this is true, then the algorithm hasnt found blocks within the acceptable range given by block_dist
    {
        
        for(j=0;j<array_num;j++)
        {
            //if the distance between them is within 3e9, to restrict number of possible calculations,  calulate the total distance between the box and photon
            #if DIMENSIONS == 2
                if ((fabs(ph_x- (*(x+j)))<block_dist) && (fabs(ph_y- (*(y+j)))<block_dist))
                {
                    
                    dist= pow(pow(ph_x- (*(x+j)), 2.0) + pow(ph_y- (*(y+j)) , 2.0),0.5);
                    //fprintf(fPtr,"Dist calculated as: %e, index: %d\n", dist, j);
                    //printf("In outer if statement, OLD: %e, %d\n", dist_min, min_index);
                    
                    if((dist<dist_min))
                    {
                        //fprintf(fPtr,"In innermost if statement, OLD: %e, %d\n", dist_min, min_index);
                        dist_min=dist; //save new minimum distance
                        min_index=j; //save index
                        //printf("New Min dist: %e, New min Index: %d, Array_Num: %d\n", dist_min, min_index, array_num);
                        
                    }
                    
                }
            #elif DIMENSIONS == 3
                if ((fabs(ph_x- (*(x+j)))<block_dist) && (fabs(ph_y- (*(y+j)))<block_dist) && (fabs(ph_z- (*(z+j)))<block_dist))
                {
                    dist= pow(pow(ph_x- (*(x+j)), 2.0) + pow(ph_y- (*(y+j)),2.0 ) + pow(ph_z- (*(z+j)) , 2.0),0.5);
                    if((dist<dist_min))
                    {
                        //printf("In innermost if statement, OLD: %e, %d\n", dist_min, min_index);
                        dist_min=dist; //save new minimum distance
                        min_index=j; //save index
                        //fprintf(fPtr,"New Min dist: %e, New min Index: %d, Array_Num: %e\n", dist_min, min_index, array_num);
                    }
                }
            #endif
        }
        block_dist*=10; //increase size of accepted distances for gris points, if dist_min==1e12 then the next time the acceptance range wil be larger
        
    }
    
    return min_index;
}

int findContainingBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z, double *szx, double *szy, int old_block_index, int find_block_switch, FILE *fPtr)
{
    int i=0, within_block_index=0;
    bool is_in_block=0; //boolean to determine if the photon is outside of a grid
    
    //can parallelize here to save time?
    for (i=0;i<array_num;i++)
    {
        
            is_in_block=checkInBlock(i,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy);
        
            if (is_in_block)
            {
                within_block_index=i;
                //change for loop index once the block is found so the code doesnt search the rest of the grids to see if the photon is within those grids
                i=array_num;
            }
        
    }
    //printf("Within Block Index:  %d\n",within_block_index);
    //if ((strcmp(DIM_SWITCH, dim_3d_str)==0) || (riken_switch==1))
    #if SIM_SWITCH == RIKEN || DIMENSIONS == 3
    {
        fprintf(fPtr, "3D switch is: %d and SIM switch is: %d\n", DIMENSIONS, SIM_SWITCH);
    }
    #endif
    
    if (is_in_block==0)
    {
        fprintf(fPtr, "Couldn't find a block that the photon is in\nx: %e y:%e\n", ph_x, ph_y);
        fflush(fPtr);
        within_block_index=-1;
    }
    
    return within_block_index;
}



int checkInBlock(int block_index, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z, double *szx, double *szy)
{
    bool is_in_block=0; //boolean to determine if the photon is outside of its previously noted block
    double x0=0, x1=0, x2=0, sz_x0=0, sz_x1=0, sz_x2=0; //coordinate and sizes of grid block, in cartesian its x,y,z in spherical its r,theta,phi
    int return_val=0;

    
        //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
        #if DIMENSIONS == 2
        {
            
            #if GEOMETRY == SPHERICAL
            {
                x0=pow(pow((*(x+block_index)),2.0)+pow((*(y+block_index)),2.0), 0.5); //radius
                x1=atan2((*(x+block_index)), (*(y+block_index))); //theta
                
                sz_x0=(*(szx+block_index));
                sz_x1=(*(szy+block_index));
                
                //pow(pow( ph_x, 2.0) + pow(ph_y, 2.0),0.5)      atan2(ph_x, ph_y)
                is_in_block= (2*fabs( ph_x - x0)- sz_x0 <= 0) && (2*fabs(ph_y - x1 ) - sz_x1 <= 0); //ph_x is ph_r for this geometry

            }
            #else
            {
                x0=(*(x+block_index));
                x1=(*(y+block_index));
                
                sz_x0=(*(szx+block_index));
                sz_x1=(*(szy+block_index));
                
                is_in_block= (2*fabs(ph_x-x0)-sz_x0 <= 0) && (2*fabs(ph_y-x1)-sz_x1 <= 0);
            }
            #endif
            
        }
        #endif
        /*
        else
        {
            if (riken_switch==1)
            {
                
                x0=pow(pow((*(x+block_index)), 2.0) + pow((*(y+block_index)),2.0 ) + pow((*(z+block_index)) , 2.0),0.5);
                x1=acos((*(z+block_index))/pow(pow((*(x+block_index)), 2.0) + pow((*(y+block_index)),2.0 ) + pow((*(z+block_index)) , 2.0),0.5));
                x2=atan2((*(y+block_index)), (*(x+block_index)));
                
                sz_x0=(*(szy+block_index));
                sz_x1=(*(szx+block_index));
                sz_x2=(*(szx+block_index));
                
                is_in_block= (fabs(pow(pow( ph_x, 2.0) + pow(ph_y, 2.0)+pow(ph_z, 2.0),0.5) - x0) <= sz_x0/2.0) &&  (fabs(acos(ph_z/pow(pow(ph_x, 2.0) + pow(ph_y,2.0 ) + pow(ph_z , 2.0),0.5)) - x1 ) <= sz_x1/2.0)  && (fabs(atan2(ph_y, ph_x) - x2 ) <= sz_x2/2.0);
                //not sure why the code was going to this line above here for spherical test
                 
             }
        }
        */
        
        if (is_in_block)
        {
            return_val=1;
        }
        else
        {
            return_val=0;
        }
    
    return return_val;
}

int findNearestPropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double hydro_domain_x, double hydro_domain_y, double epsilon_b, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *all_time_steps, int *sorted_indexes, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr)
{
    int i=0, min_index=0, ph_block_index=0;
    double ph_x=0, ph_y=0, ph_phi=0, ph_z=0, ph_r=0, ph_theta=0;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates

    double ph_v_norm=0, fl_v_norm=0, synch_x_sect=0;
    double n_cosangle=0, n_dens_lab_tmp=0,n_vx_tmp=0, n_vy_tmp=0, n_vz_tmp=0, n_temp_tmp=0 ;
    double rnd_tracker=0, n_dens_min=0, n_vx_min=0, n_vy_min=0, n_vz_min=0, n_temp_min=0;
    int num_thread=omp_get_num_threads();
    bool is_in_block=0; //boolean to determine if the photon is outside of its previously noted block
    
    int index=0, num_photons_find_new_element=0;
    double mfp=0,min_mfp=0, beta=0;
    double el_p[4];
    double ph_p_comv[4], ph_p[4], fluid_beta[3];

    //initialize gsl random number generator fo each thread
    
        const gsl_rng_type *rng_t;
        gsl_rng **rng;
        gsl_rng_env_setup();
        rng_t = gsl_rng_ranlxs0;

        rng = (gsl_rng **) malloc((num_thread ) * sizeof(gsl_rng *)); 
        rng[0]=rand;

            //#pragma omp parallel for num_threads(nt)
        for(i=1;i<num_thread;i++)
        {
            rng[i] = gsl_rng_alloc (rng_t);
            gsl_rng_set(rng[i],gsl_rng_get(rand));
        }
       
    //go through each photon and find the blocks around it and then get the distances to all of those blocks and choose the one thats the shortest distance away
    //can optimize here, exchange the for loops and change condition to compare to each of the photons is the radius of the block is .95 (or 1.05) times the min (max) photon radius
    //or just parallelize this part here
    
    min_mfp=1e12;
    #pragma omp parallel for num_threads(num_thread) firstprivate( is_in_block, ph_block_index, ph_x, ph_y, ph_z, ph_phi, ph_r, min_index, n_dens_lab_tmp,n_vx_tmp, n_vy_tmp, n_vz_tmp, n_temp_tmp, fl_v_x, fl_v_y, fl_v_z, fl_v_norm, ph_v_norm, n_cosangle, mfp, beta, rnd_tracker, ph_p_comv, el_p, ph_p, fluid_beta) private(i) shared(min_mfp ) reduction(+:num_photons_find_new_element)
    for (i=0;i<num_ph; i++)
    {
        //fprintf(fPtr, "%d, %d,%e\n", i, ((ph+i)->nearest_block_index), ((ph+i)->weight));
        //fflush(fPtr);
        
        if (find_nearest_block_switch==0)
        {
            ph_block_index=(ph+i)->nearest_block_index; //if starting a new frame the number of indexes can change and cause a seg fault here
        }
        else
        {
            ph_block_index=0; // therefore if starting a new frame set index=0 to avoid this issue
        }
        
        //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
        #if DIMENSIONS == 2
        {
            #if GEOMETRY == SPHERICAL
                ph_x=pow(pow(((ph+i)->r0),2.0)+pow(((ph+i)->r1),2.0), 0.5); //convert back to 2d spherical coordinate
                ph_y=((ph+i)->r2);
                ph_r=pow(ph_x*ph_x + ph_y*ph_y, 0.5);
                ph_theta=acos(ph_y/ph_r); //this is actually theta in this context
                ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
            #elif GEOMETRY == CARTESIAN
                ph_x=pow(pow(((ph+i)->r0),2.0)+pow(((ph+i)->r1),2.0), 0.5); //convert back to FLASH x coordinate (2d cartesian hydro coordinates)
                ph_y=((ph+i)->r2);
                ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
                ph_r=pow(ph_x*ph_x + ph_y*ph_y, 0.5);
            #endif
            
        }
        #else
        {
            ph_x=((ph+i)->r0);
            ph_y=((ph+i)->r1);
            ph_z=((ph+i)->r2);
            ph_r=pow(ph_x*ph_x + ph_y*ph_y+ph_z*ph_z, 0.5);
        }
        #endif
        //printf("ph_x:%e, ph_y:%e\n", ph_x, ph_y);
        
        //if the location of the photon is less than the domain of the hydro simulation then do all of this, otherwise assing huge mfp value so no scattering occurs and the next frame is loaded
        // absorbed photons have ph_block_index=-1, therefore if this value is not less than 0, calulate the mfp properly but doesnt work when go to new frame and find new indexes (will change b/c will get rid of these photons when printing)
        //alternatively make decision based on 0 weight
        if (((ph_y<hydro_domain_y) && (ph_x<hydro_domain_x)) && ((ph+i)->nearest_block_index != -1) ) //can use sorted index to see which photons have been absorbed efficiently before printing and get the indexes
        {
            #if GEOMETRY == SPHERICAL
                is_in_block=checkInBlock(ph_block_index,  ph_r,  ph_theta,  ph_z,  x,   y, z,  szx,  szy);
            #elif GEOMETRY == CARTESIAN
                is_in_block=checkInBlock(ph_block_index,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy);
            #endif
        
            if (find_nearest_block_switch==0 && is_in_block)
            {
                //keep the saved grid index
                min_index=ph_block_index;
            }
            else
            {
                //find the new index of the block closest to the photon
                //min_index=findNearestBlock(array_num,  ph_x,  ph_y,  ph_z,  x,   y,  z); //stop doing this one b/c nearest grid could be one that the photon isnt actually in due to adaptive mesh
            
                //find the new index of the block that the photon is actually in
                #if DIMENSIONS == 2
                {
                    #if GEOMETRY == SPHERICAL
                        min_index=findContainingBlock(array_num,  ph_r,  ph_theta,  ph_z,  x,   y, z,  szx,  szy, ph_block_index, find_nearest_block_switch, fPtr);
                    #elif GEOMETRY == CARTESIAN
                        min_index=findContainingBlock(array_num,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy, ph_block_index, find_nearest_block_switch, fPtr);
                    #endif
                    }
                #endif
                if (min_index != -1)
                {
                    (ph+i)->nearest_block_index=min_index; //save the index if min_index != -1
                    
                    //also recalculate the photons' comoving frequency in this new fluid element
                    ph_p[0]=((ph+i)->p0);
                    ph_p[1]=((ph+i)->p1);
                    ph_p[2]=((ph+i)->p2);
                    ph_p[3]=((ph+i)->p3);
                    
                    //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                    #if DIMENSIONS == 2
                    {
                        fluid_beta[0]=(*(velx+min_index))*cos(ph_phi);
                        fluid_beta[1]=(*(velx+min_index))*sin(ph_phi);
                        fluid_beta[2]=(*(vely+min_index));
                    }
                    #else
                    {
                        fluid_beta[0]=(*(velx+min_index));
                        fluid_beta[1]=(*(vely+min_index));
                        fluid_beta[2]=(*(velz+min_index));
                    }
                    #endif
                    
                    lorentzBoost(&fluid_beta, &ph_p, &ph_p_comv, 'p', fPtr);
                    
                    ((ph+i)->comv_p0)=ph_p_comv[0];
                    ((ph+i)->comv_p1)=ph_p_comv[1];
                    ((ph+i)->comv_p2)=ph_p_comv[2];
                    ((ph+i)->comv_p3)=ph_p_comv[3];
                    
                    num_photons_find_new_element+=1;
                }
                else
                {
                	fprintf(fPtr, "Photon number %d FLASH index not found, making sure it doesnt scatter.\n", i);
                }
            
            }
            
            //if min_index!= -1 (know which fluid element photon is in) do all this stuff, otherwise make sure photon doesnt scatter
            if (min_index != -1)
            {
                //fprintf(fPtr,"Min Index: %d\n", min_index);
        
                //save values
                (n_dens_lab_tmp)= (*(dens_lab+min_index));
                (n_vx_tmp)= (*(velx+min_index));
                (n_vy_tmp)= (*(vely+min_index));
                (n_temp_tmp)= (*(temp+min_index));
                

                //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
                #if DIMENSIONS == 3
                {
                    (n_vz_tmp)= (*(velz+min_index));
                }
                #endif
        
                //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                #if DIMENSIONS == 2
                {
                    fl_v_x=(*(velx+min_index))*cos(ph_phi);
                    fl_v_y=(*(velx+min_index))*sin(ph_phi);
                    fl_v_z=(*(vely+min_index));
                }
                #else
                {
                    fl_v_x=(*(velx+min_index));
                    fl_v_y=(*(vely+min_index));
                    fl_v_z=(*(velz+min_index));
                }
                #endif
                
                fl_v_norm=pow(pow(fl_v_x, 2.0)+pow(fl_v_y, 2.0)+pow(fl_v_z, 2.0), 0.5);
                ph_v_norm=pow(pow(((ph+i)->p1), 2.0)+pow(((ph+i)->p2), 2.0)+pow(((ph+i)->p3), 2.0), 0.5);
        
                //(*(n_cosangle+i))=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //find cosine of the angle between the photon and the fluid velocities via a dot product
                (n_cosangle)=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //make 1 for cylindrical otherwise its undefined
        
                //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                #if DIMENSIONS == 2
                {
                    beta=pow((n_vx_tmp*n_vx_tmp)+(n_vy_tmp*n_vy_tmp),0.5);
                }
                #else
                {
                    beta=pow((pow((n_vx_tmp),2)+pow((n_vy_tmp),2)+pow((n_vz_tmp),2)),0.5);
                }
                #endif

                *(ph_p+0)=((ph+i)->p0);
                *(ph_p+1)=((ph+i)->p1);
                *(ph_p+2)=((ph+i)->p2);
                *(ph_p+3)=((ph+i)->p3);
                
                //ph_p_comv[0]=((ph+i)->comv_p0);
                //ph_p_comv[1]=((ph+i)->comv_p1);
                //ph_p_comv[2]=((ph+i)->comv_p2);
                //ph_p_comv[3]=((ph+i)->comv_p3);
                
                //printf("ph: p0 %e p1 %e p2 %e p3 %e\n",  *(ph_p_comv+0), *(ph_p_comv+1), *(ph_p_comv+2), *(ph_p_comv+3));

                
                //singleElectron(&el_p[0], n_temp_tmp, &ph_p_comv[0], rng[omp_get_thread_num()], fPtr); //get random electron
                //printf("after singleElectron n_temp_tmp %e from ptr %e n_dens_tmp %e from ptr %e\n", n_temp_tmp, (*(temp+min_index)), n_dens_tmp, (*(dens+min_index)));
                
                //printf("Chosen el: p0 %e p1 %e p2 %e p3 %e\nph: p0 %e p1 %e p2 %e p3 %e\n", *(el_p+0), *(el_p+1), *(el_p+2), *(el_p+3), *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3));
                
                //synch_x_sect=synCrossSection(n_dens_tmp/M_P, n_temp_tmp, ph_p_comv[0]*C_LIGHT/PL_CONST, sqrt((el_p[0]*el_p[0]/(M_EL*M_EL*C_LIGHT*C_LIGHT))-1), epsilon_b);
                //printf("i: %d flash_array_idx %d synch_x_sect %e freq %e temp %e el_dens %e\n", i, min_index, synch_x_sect, *(ph_p+0)*C_LIGHT/PL_CONST, n_temp_tmp, n_dens_tmp/M_P);
                
                //if (synch_x_sect==0)
                //{
                //*(will_scatter+i)=1; //this photon will scatter b/c probability of absorption=0
                //}
                /*
                else
                {
                    if (gsl_rng_uniform_pos(rng[omp_get_thread_num()])>(THOM_X_SECT/(THOM_X_SECT+synch_x_sect)))
                    {
                        //this photon will be absorbed
                        *(will_scatter+i)=0;
                    }
                    else
                    {
                        *(will_scatter+i)=1;
                    }
                    
                } photons can onlt scatter now
                */
                
                //put this in to double check that random number is between 0 and 1 (exclusive) because there was a problem with this for parallel case
                rnd_tracker=0;
        
                rnd_tracker=gsl_rng_uniform_pos(rng[omp_get_thread_num()]);
                //printf("Rnd_tracker: %e Thread number %d \n",rnd_tracker, omp_get_thread_num() );
        
                //mfp=(-1)*log(rnd_tracker)*(M_P/((n_dens_tmp))/(THOM_X_SECT)); ///(1.0-beta*((n_cosangle)))) ; //calulate the mfp and then multiply it by the ln of a random number to simulate distribution of mean free paths DO EVERYTHING IN COMOV FRAME NOW
                mfp=(-1)*(M_P/((n_dens_lab_tmp))/THOM_X_SECT/(1.0-beta*((n_cosangle))))*log(rnd_tracker) ;
                
                //if (mfp/C_LIGHT < 1e-100)
                //{
                //    fprintf("Photon %d has a mfp of %d\n", i, mfp);
                //    exit(0);
                //}
                
            }
            else
            {
                mfp=min_mfp;
            }
        }
         else
        {
            mfp=min_mfp;
            //fprintf(fPtr,"Photon %d In ELSE\n", i);
            //exit(0);
        }
        
        *(all_time_steps+i)=mfp/C_LIGHT;
    }
    //exit(0);
    //free rand number generator
    for (i=1;i<num_thread;i++)
    {
        gsl_rng_free(rng[i]);
    }
    free(rng);
        
    //printf("HERE\n");
    for (i=0;i<num_ph;i++)
    {
        *(sorted_indexes+i)= i; //save  indexes to array to use in qsort
    }
    
    //printf("before QSORT\n");
    #if (defined _GNU_SOURCE || defined __GNU__ || defined __linux__)
        qsort_r(sorted_indexes, num_ph, sizeof (int),  compare2, all_time_steps);
    #elif (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || defined __FREEBSD__ || defined __BSD__ || defined OpenBSD3_1 || defined OpenBSD3_9)
        qsort_r(sorted_indexes, num_ph, sizeof (int), all_time_steps, compare);
    #else
        #error Cannot detect operating system
    #endif
    
    //for (i=0;i<num_ph;i++)
    //{
    //    fprintf(fPtr, "Qsort: %d GSL: %d\n", *(sorted_indexes_2+i), *(sorted_indexes+i));
    //}
    //exit(0);
    
    //print number of times we had to refind the index of the elemtn photons were located in
    if (find_nearest_block_switch!=0)
    {
        num_photons_find_new_element=0; //force this to be 0 since we forced MCRaT to find the indexes for all the photons here
    }
    
    //fprintf(fPtr, "MCRat had to refind where %d photons were located in the grid\n", num_photons_find_new_element);
    //(*time_step)=*(all_time_steps+(*(sorted_indexes+0))); //dont need to return index b/c photonEvent doesnt use this, but mcrat.c uses this info
    //index= *(sorted_indexes+0);//first element of sorted array
    //free(el_p);free(ph_p_comv);
    return num_photons_find_new_element;
    
}

int compare (void *ar, const void *a, const void *b)
{
    //from https://phoxis.org/2012/07/12/get-sorted-index-orderting-of-an-array/
  int aa = *(int *) a;
  int bb = *(int *) b;
  double *arr=NULL;
  arr=ar;
  
  //printf("%d, %d\n", aa, bb);
  //printf("%e, %e\n", arr[aa] , arr[bb]);
  //return (aa - bb);
  /*
 if (arr[aa] < arr[bb])
    return -1; 
  if (arr[aa] == arr[bb])
    return 0;
  if (arr[aa] > arr[bb])
    return 1;
    */
    return ((arr[aa] > arr[bb]) - (arr[aa] < arr[bb]));
}

int compare2 ( const void *a, const void *b, void *ar)
{
    //have 2 compare funcions b/c of changes in qsort_r between BSD and GNU
    //from https://phoxis.org/2012/07/12/get-sorted-index-orderting-of-an-array/
  int aa = *(int *) a;
  int bb = *(int *) b;
  double *arr=NULL;
  arr=ar;
  
  //printf("%d, %d\n", aa, bb);
  //printf("%e, %e\n", arr[aa] , arr[bb]);
  //return (aa - bb);
  /*
 if (arr[aa] < arr[bb])
    return -1; 
  if (arr[aa] == arr[bb])
    return 0;
  if (arr[aa] > arr[bb])
    return 1;
    */
    return ((arr[aa] > arr[bb]) - (arr[aa] < arr[bb]));
}

int interpolatePropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *n_dens_lab, double *n_vx, double *n_vy, double *n_vz, double *n_temp, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr)
{
    /*
     * THIS FUNCTION IS WRITTEN JUST FOR 2D SIMS AS OF NOW, not used
    */
    int i=0, j=0, min_index=0, ph_block_index=0;
    int left_block_index=0, right_block_index=0, bottom_block_index=0, top_block_index=0, all_adjacent_block_indexes[4];
    double ph_x=0, ph_y=0, ph_phi=0, ph_z=0, dist=0, left_dist_min=0, right_dist_min=0, top_dist_min=0, bottom_dist_min=0, dv=0, v=0;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates
    double r=0, theta=0;

    double ph_v_norm=0, fl_v_norm=0;
    double n_cosangle=0, n_dens_lab_tmp=0,n_vx_tmp=0, n_vy_tmp=0, n_vz_tmp=0, n_temp_tmp=0;
    double rnd_tracker=0, n_dens_lab_min=0, n_vx_min=0, n_vy_min=0, n_vz_min=0, n_temp_min=0;
    int num_thread=2;//omp_get_max_threads();
    bool is_in_block=0; //boolean to determine if the photon is outside of its previously noted block
    
    int index=0;
    double mfp=0,min_mfp=0, beta=0;
        
        
    //initialize gsl random number generator fo each thread
    
        const gsl_rng_type *rng_t;
        gsl_rng **rng;
        gsl_rng_env_setup();
        rng_t = gsl_rng_ranlxs0;

        rng = (gsl_rng **) malloc((num_thread ) * sizeof(gsl_rng *)); 
        rng[0]=rand;

            //#pragma omp parallel for num_threads(nt)
        for(i=1;i<num_thread;i++)
        {
            rng[i] = gsl_rng_alloc (rng_t);
            gsl_rng_set(rng[i],gsl_rng_get(rand));
        }
       
    //go through each photon and find the blocks around it and then get the distances to all of those blocks and choose the one thats the shortest distance away
    //can optimize here, exchange the for loops and change condition to compare to each of the photons is the radius of the block is .95 (or 1.05) times the min (max) photon radius
    //or just parallelize this part here
    
    min_mfp=1e12;
    #pragma omp parallel for num_threads(num_thread) firstprivate( r, theta,dv, v, all_adjacent_block_indexes, j, left_block_index, right_block_index, top_block_index, bottom_block_index, is_in_block, ph_block_index, ph_x, ph_y, ph_z, ph_phi, min_index, n_dens_lab_tmp,n_vx_tmp, n_vy_tmp, n_vz_tmp, n_temp_tmp, fl_v_x, fl_v_y, fl_v_z, fl_v_norm, ph_v_norm, n_cosangle, mfp, beta, rnd_tracker) private(i) shared(min_mfp )
    for (i=0;i<num_ph; i++)
    {
        //printf("%d, %e,%e\n", i, ((ph+i)->r0), ((ph+i)->r1));
        if (find_nearest_block_switch==0)
        {
            ph_block_index=(ph+i)->nearest_block_index; //if starting a new frame the number of indexes can change and cause a seg fault
        }
        else
        {
            ph_block_index=0; //if starting a new frame set index=0 to avoid this issue
        }
        
        //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
        #if DIMENSIONS == 2
        {
            ph_x=pow(pow(((ph+i)->r0),2.0)+pow(((ph+i)->r1),2.0), 0.5); //convert back to FLASH x coordinate
            ph_y=((ph+i)->r2);
            ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
            
        }
        #else
        {
            ph_x=((ph+i)->r0);
            ph_y=((ph+i)->r1);
            ph_z=((ph+i)->r2);
            
        }
        #endif
        //printf("ph_x:%e, ph_y:%e\n", ph_x, ph_y);
        
        is_in_block=checkInBlock(ph_block_index,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy);
        
        if (find_nearest_block_switch==0 && is_in_block)
        {
            //keep the saved grid index
            min_index=ph_block_index;
        }
        else
        {
            //find the new index of the block closest to the photon
            //min_index=findNearestBlock(array_num,  ph_x,  ph_y,  ph_z,  x,   y,  z); //stop doing this one b/c nearest grid could be one that the photon isnt actually in due to adaptive mesh
            
            //find the new index of the block that the photon is actually in
            min_index=findContainingBlock(array_num,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy, ph_block_index, find_nearest_block_switch, fPtr);
            
            (ph+i)->nearest_block_index=min_index; //save the index
            
        }
        
        //look for the blocks surounding the block of interest and order them by the 
        left_dist_min=1e15;//set dist to impossible value to make sure at least first distance calulated is saved
        right_dist_min=1e15;
        top_dist_min=1e15;
        bottom_dist_min=1e15;
        for (j=0;j<array_num;j++)
        {
            //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
            #if DIMENSIONS == 2
            {
                dist= pow(pow((*(x+min_index))- (*(x+j)), 2.0) + pow((*(y+min_index))- (*(y+j)) , 2.0),0.5);
            }
            #else
            {
                dist= pow(pow((*(x+min_index))- (*(x+j)), 2.0) + pow((*(y+min_index))- (*(y+j)),2.0 ) + pow((*(z+min_index))- (*(z+j)) , 2.0),0.5);
            }
            #endif
            
            if ((*(x+j))<(*(x+min_index)) && (dist < left_dist_min) )
            {
                left_block_index=j;
                left_dist_min=dist;
            }
            else if ((*(x+j))>(*(x+min_index)) && (dist < right_dist_min))
            {
                right_block_index=j;
                right_dist_min=dist;
            }
            
            if ((*(y+j))<(*(y+min_index)) && (dist < bottom_dist_min) )
            {
                bottom_block_index=j;
                bottom_dist_min=dist;
            }
            else if ((*(y+j))>(*(y+min_index)) && (dist < top_dist_min) )
            {
                top_block_index=j;
                top_dist_min=dist;
            }
        
        }
        all_adjacent_block_indexes[0]=left_block_index;
        all_adjacent_block_indexes[1]=right_block_index;
        all_adjacent_block_indexes[2]=bottom_block_index;
        all_adjacent_block_indexes[3]=top_block_index;       
        
        //do a weighted average of the 4 nearest grids based on volume
        v=0;
        (n_dens_lab_tmp)=0;
        (n_vx_tmp)= 0;
        (n_vy_tmp)= 0;
        (n_temp_tmp)= 0;
        (n_vz_tmp)= 0;
            
        for (j=0;j<4;j++)
        {
            
            #if SIM_SWITCH == RIKEN
            {
                r=pow(pow((*(x+all_adjacent_block_indexes[j])),2.0)+pow((*(y+all_adjacent_block_indexes[j])),2.0), 0.5);
                theta=atan2((*(x+all_adjacent_block_indexes[j])), (*(y+all_adjacent_block_indexes[j])));
                dv=2.0*M_PI*pow(r,2)*sin(theta)*(*(szx+all_adjacent_block_indexes[j]))*(*(szy+all_adjacent_block_indexes[j])) ;
            }
            #else
            {
                //using FLASH
                dv=2.0*M_PI*(*(x+all_adjacent_block_indexes[j]))*pow(*(szx+all_adjacent_block_indexes[j]),2.0)  ;

            }
            #endif
            
            v+=dv;
            
            //save values
            (n_dens_lab_tmp)+= (*(dens_lab+all_adjacent_block_indexes[j]))*dv;
            (n_vx_tmp)+= (*(velx+all_adjacent_block_indexes[j]))*dv;
            (n_vy_tmp)+= (*(vely+all_adjacent_block_indexes[j]))*dv;
            (n_temp_tmp)+= (*(temp+all_adjacent_block_indexes[j]))*dv;
            
            //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
            #if DIMENSIONS == 3
            {
                (n_vz_tmp)+= (*(velz+all_adjacent_block_indexes[j]))*dv;
            }
            #endif
            
        }
        

         //fprintf(fPtr,"Outside\n");
        
        //save values
        (n_dens_lab_tmp)/= v;
        (n_vx_tmp)/= v;
        (n_vy_tmp)/= v;
        (n_temp_tmp)/= v;
        //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
        #if DIMENSIONS == 3
        {
            (n_vz_tmp)/= v;
        }
        #endif
        
        //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
        #if DIMENSIONS == 2
        {
            fl_v_x=n_vx_tmp*cos(ph_phi);
            fl_v_y=n_vx_tmp*sin(ph_phi);
            fl_v_z=n_vy_tmp;
        }
        #else
        {
            fl_v_x=n_vx_tmp;
            fl_v_y=n_vy_tmp;
            fl_v_z=n_vz_tmp;
        }
        #endif
        
        fl_v_norm=pow(pow(fl_v_x, 2.0)+pow(fl_v_y, 2.0)+pow(fl_v_z, 2.0), 0.5);
        ph_v_norm=pow(pow(((ph+i)->p1), 2.0)+pow(((ph+i)->p2), 2.0)+pow(((ph+i)->p3), 2.0), 0.5);
        
        //(*(n_cosangle+i))=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //find cosine of the angle between the photon and the fluid velocities via a dot product
        (n_cosangle)=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //make 1 for cylindrical otherwise its undefined
        
        //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
        #if DIMENSIONS == 2
        {
            beta=pow((pow((n_vx_tmp),2)+pow((n_vy_tmp),2)),0.5);
        }
        #else
        {
            beta=pow((pow((n_vx_tmp),2)+pow((n_vy_tmp),2)+pow((n_vz_tmp),2)),0.5);
        }
        #endif
        //put this in to double check that random number is between 0 and 1 (exclusive) because there was a problem with this for parallel case
        rnd_tracker=0;
        
        rnd_tracker=gsl_rng_uniform_pos(rng[omp_get_thread_num()]);
        
        mfp=(-1)*(M_P/((n_dens_lab_tmp))/THOM_X_SECT/(1.0-beta*((n_cosangle))))*log(rnd_tracker) ; //calulate the mfp and then multiply it by the ln of a random number to simulate distribution of mean free paths 
        
        
        #pragma omp critical 
        if ( mfp<min_mfp)
        {
            min_mfp=mfp;
            n_dens_lab_min= n_dens_lab_tmp;
            n_vx_min= n_vx_tmp;
            n_vy_min= n_vy_tmp;
            //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
            #if DIMENSIONS == 3
            {
                n_vz_min= n_vz_tmp;
            }
            #endif
            
            n_temp_min= n_temp_tmp;
            index=i;
            //fprintf(fPtr, "Thread is %d. new min: %e for photon %d with block properties: %e, %e, %e Located at: %e, %e, Dist: %e\n", omp_get_thread_num(), mfp, index, n_vx_tmp, n_vy_tmp, n_temp_tmp, *(x+min_index), *(y+min_index), dist_min);
            //fflush(fPtr);
            #pragma omp flush(min_mfp)
        }

        
    }
    
    //free rand number generator
    for (i=1;i<num_thread;i++)
    {
        gsl_rng_free(rng[i]);
    }
    free(rng);
    
    *(n_dens_lab)= n_dens_lab_min;
    *(n_vx)= n_vx_min;
    *(n_vy)= n_vy_min;
    //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
    #if DIMENSIONS == 3
    {
        *(n_vz)= n_vz_min;
    }
    #endif
    
    *(n_temp)= n_temp_min;
    (*time_step)=min_mfp/C_LIGHT;
    return index;
    
}


void updatePhotonPosition(struct photon *ph, int num_ph, double t, FILE *fPtr)
{
    //move photons by speed of light
 
    int i=0, num_thread=omp_get_num_threads();
    double old_position=0, new_position=0, divide_p0=0;
    
    
    #pragma omp parallel for num_threads(num_thread) firstprivate(old_position, new_position, divide_p0)
    for (i=0;i<num_ph;i++)
    {
        if (((ph+i)->type != 's') && ((ph+i)->weight != 0))
        {
            old_position= pow(  pow((ph+i)->r0,2)+pow((ph+i)->r1,2)+pow((ph+i)->r2,2), 0.5 ); //uncommented checks since they were not necessary anymore
            
            divide_p0=1.0/((ph+i)->p0);
            
            ((ph+i)->r0)+=((ph+i)->p1)*divide_p0*C_LIGHT*t; //update x position
            
            ((ph+i)->r1)+=((ph+i)->p2)*divide_p0*C_LIGHT*t;//update y
            
            ((ph+i)->r2)+=((ph+i)->p3)*divide_p0*C_LIGHT*t;//update z
            
            new_position= pow(  pow((ph+i)->r0,2)+pow((ph+i)->r1,2)+pow((ph+i)->r2,2), 0.5 );
            
            //if ((new_position-old_position)/t > C_LIGHT)
            //{
            //    fprintf(fPtr, "PHOTON NUMBER %d IS SUPERLUMINAL. ITS SPEED IS %e c.\n", i, ((new_position-old_position)/t)/C_LIGHT);
            //}
            
            if ( (ph+i)->s0 != 1)
            {
            //	fprintf(fPtr, "PHOTON NUMBER %d DOES NOT HAVE I=1. Instead it is: %e\n", i, (ph+i)->s0);
            }
            
            //printf("In update  function: %e, %e, %e, %e, %e, %e, %e\n",((ph+i)->r0), ((ph+i)->r1), ((ph+i)->r2), t, ((ph+i)->p1)/((ph+i)->p0), ((ph+i)->p2)/((ph+i)->p0), ((ph+i)->p3)/((ph+i)->p0) );
        }
    }
        
    //printf("In update  function: %e, %e, %e, %e\n",t, ((ph)->p1)/((ph)->p0), ((ph)->p2)/((ph)->p0), ((ph)->p3)/((ph)->p0) );    
    
}

void mullerMatrixRotation(double theta, double *s, FILE *fPtr)
{
    //makes a CCW rotation od the stokes parameters when the photon velocity vector is pointed towards the observer, follows Lundman
    gsl_matrix *M= gsl_matrix_calloc (4, 4); //create matrix thats 4x4 to do rotation as defined in McMaster 1961 (has it to rotate CW in that paper)
    gsl_vector *result= gsl_vector_alloc(4);
    gsl_vector_view stokes;
    
    stokes=gsl_vector_view_array(s, 4);
    //fprintf(fPtr, "sokes parameter before= %e %e %e %e\n", gsl_vector_get(&stokes.vector, 0), gsl_vector_get(&stokes.vector, 1), gsl_vector_get(&stokes.vector, 2), gsl_vector_get(&stokes.vector, 3));
    
    gsl_matrix_set(M, 0,0,1);
    gsl_matrix_set(M, 3,3,1);
    gsl_matrix_set(M, 1,1,cos(2*theta));
    gsl_matrix_set(M, 2,2,cos(2*theta));
    gsl_matrix_set(M, 1,2,-1*sin(2*theta));
    gsl_matrix_set(M, 2,1,sin(2*theta));
    gsl_blas_dgemv(CblasNoTrans, 1, M, &stokes.vector, 0, result); //Ms=s
    
    //fprintf(fPtr, "stokes parameter after= %e %e %e %e\n\n", gsl_vector_get(result, 0), gsl_vector_get(result, 1), gsl_vector_get(result, 2), gsl_vector_get(result, 3));
    
    //save back to the original stokes vector
    *(s+0)=gsl_vector_get(result, 0);
    *(s+1)=gsl_vector_get(result, 1);
    *(s+2)=gsl_vector_get(result, 2);
    *(s+3)=gsl_vector_get(result, 3);
    
    gsl_vector_free(result);
    gsl_matrix_free (M);
    
}

void findXY(double *v_ph, double *vector, double *x, double *y)
{
    //finds the stokes plane coordinate x,y axis for the photon velocity with respect to some reference vector
    //assumes that pointers point to array of 3 doubles in length
    double norm=0;
    
    *(y+0)= ((*(v_ph+1))*(*(vector+2))-(*(v_ph+2))*(*(vector+1)));
    *(y+1)= -1*((*(v_ph+0))*(*(vector+2))-(*(v_ph+2))*(*(vector+0)));
    *(y+2)= ((*(v_ph+0))*(*(vector+1))-(*(v_ph+1))*(*(vector+0))); // vector X v_ph
    
    norm=1.0/sqrt( (*(y+0))*(*(y+0)) + (*(y+1))*(*(y+1)) + (*(y+2))*(*(y+2)));
    *(y+0) *= norm;
    *(y+1) *= norm;
    *(y+2) *= norm;
    
    *(x+0)= (*(y+1))*(*(v_ph+2))-(*(y+2))*(*(v_ph+1));
    *(x+1)= -1*((*(y+0))*(*(v_ph+2))-(*(y+2))*(*(v_ph+0)));
    *(x+2)= (*(y+0))*(*(v_ph+1))-(*(y+1))*(*(v_ph+0));
    
    norm=1.0/sqrt( (*(x+0))*(*(x+0)) + (*(x+1))*(*(x+1)) + (*(x+2))*(*(x+2)));
    *(x+0) *= norm;
    *(x+1) *= norm;
    *(x+2) *= norm;
    
}

double findPhi(double *x_old, double *y_old, double *x_new, double *y_new)
{
    //find the angle to rotate the stokes vector to transform from one set of stokes coordinates to another
    //this is given by Lundman
    gsl_vector_view y=gsl_vector_view_array(y_old, 3);
    gsl_vector_view x=gsl_vector_view_array(x_old, 3);
    gsl_vector_view y_prime=gsl_vector_view_array(y_new, 3);
    gsl_vector_view x_prime=gsl_vector_view_array(x_new, 3);
    double factor=0, dot_prod_result=0;
    
    gsl_blas_ddot(&x.vector, &y_prime.vector, &dot_prod_result);
    
    if (dot_prod_result>0)
    {
        factor=1;
    }
    else if (dot_prod_result<0)
    {
        factor=-1;
    }
    else
    {
        factor=0;
    }
    
    gsl_blas_ddot(&y.vector, &y_prime.vector, &dot_prod_result);
    
    return -1*factor*acos(dot_prod_result);
}

void stokesRotation(double *v, double *v_ph, double *v_ph_boosted, double *s, FILE *fPtr)
{
    //takes 3 velocities of the initial photon, v_ph, the boosted photon, v_ph_boosted. and the boost vector, v
    double z_hat[3]={0,0,1}; //z to calulate stokes
    double x[3]={0,0,0}, y[3]={0,0,0}, x_new[3]={0,0,0}, y_new[3]={0,0,0};//initalize arrays to hold stokes coordinate system
    double phi=0;
    
    //if (i==0)
    {
    //find stokes coordinate sys in orig frame with respect to z axis
    findXY(v_ph, &z_hat, &x, &y);
    }
    
    //find stokes coordinate sys in orig frame with respect to boost vector
    findXY(v_ph, v, &x_new, &y_new);
    
    phi=findPhi(x, y, x_new, y_new);//now find rotation between the two coordinate systems
    
    //rotate the stokes vector now to put it in the coordinate system fo the boosted photon and the boost evctor
    mullerMatrixRotation(phi, s, fPtr);
    
    //find the new coordinates of the rotated stokes vector with the boosted photon and the boost vector
    findXY(v_ph_boosted, v, &x, &y);
    
    //find stokes coordinate sys in orig frame with respect to z axis
    findXY(v_ph_boosted, &z_hat, &x_new, &y_new);
    
    phi=findPhi(x, y, x_new, y_new);//now find rotation between the two coordinate systems
    
    //do the rotation of the stokes vector to put it in the coordinate system of the boosted photon and the z axis
    mullerMatrixRotation(phi, s, fPtr);
    
}


double photonEvent(struct photon *ph, int num_ph, double dt_max, double *all_time_steps, int *sorted_indexes, double *all_flash_vx, double *all_flash_vy, double *all_flash_vz, double *all_fluid_temp, int *scattered_ph_index, int *frame_scatt_cnt, int *frame_abs_cnt, gsl_rng * rand, FILE *fPtr)
{
    //function to perform single photon scattering
    int  i=0, index=0, ph_index=0, event_did_occur=0; //variable event_did_occur is to keep track of wether a scattering or absorption actually occured or not,
    double scatt_time=0, old_scatt_time=0; //keep track of new time to scatter vs old time to scatter to know how much to incrementally propagate the photons if necessary
    double phi=0, theta=0; //phi and theta for the 4 momentum 
    double ph_phi=0, flash_vx=0, flash_vy=0, flash_vz=0, fluid_temp=0;    
    double *ph_p=malloc(4*sizeof(double)); //pointer to hold only photon 4 momentum @ start
    double *el_p_comov=malloc(4*sizeof(double));//pointer to hold the electron 4 momenta in comoving frame
    double *ph_p_comov=malloc(4*sizeof(double));//pointer to hold the comoving photon 4 momenta
    double *fluid_beta=malloc(3*sizeof(double));//pointer to hold fluid velocity vector
    double *negative_fluid_beta=malloc(3*sizeof(double));//pointer to hold negative fluid velocity vector
    double *s=malloc(4*sizeof(double)); //vector to hold the stokes parameters for a given photon
    
    i=0;
    old_scatt_time=0;
    event_did_occur=0;
    //fprintf(fPtr,"In this function Num_ph %d\n", num_ph);
    //fflush(fPtr);
        
    while (i<num_ph && event_did_occur==0 )
    {
        ph_index=(*(sorted_indexes+i));
        
        scatt_time= *(all_time_steps+ph_index); //get the time until the photon scatters
        
        //IF THE TIME IS GREATER THAN dt_max dont let the photons positions be updated
        if (scatt_time<dt_max)
        {
            updatePhotonPosition(ph, num_ph, scatt_time-old_scatt_time, fPtr);
        
            //fprintf(fPtr,"i: %d, Photon: %d, Delta t=%e\n", i, ph_index, scatt_time-old_scatt_time);
            //fflush(fPtr);
            
            //if the photon should scatter then do so, will_scatter==1
            //if (*(will_scatter+ph_index) != 0 ) ont need b/c all photns are able to scatter and none can be explicitly absorbed
            //{
            
                //WHAT IF THE PHOTON MOVES TO A NEW BLOCK BETWEEN WHEN WE CALC MFP AND MOVE IT TO DO THE SCATTERING????
                //it mostly happens at low optical depth, near the photosphere so we would have a large mfp anyways so we probably wouldn't be in this function in that case
                index=(ph+ph_index)->nearest_block_index; //the sorted_indexes gives index of photon with smallest time to potentially scatter then extract the index of the block closest to that photon
        
                flash_vx=*(all_flash_vx+  index);
                flash_vy=*(all_flash_vy+  index);
                fluid_temp=*(all_fluid_temp+  index);
                //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
                #if DIMENSIONS == 3
                {
                    flash_vz=*(all_flash_vz+  index);
                }
                #endif
        
                ph_phi=atan2(((ph+ph_index)->r1), (((ph+ph_index)->r0)));
                /*
                fprintf(fPtr,"ph_phi=%e\n", ph_phi);
                fflush(fPtr);
                */

                //convert flash coordinated into MCRaT coordinates
                //printf("Getting fluid_beta\n");
        
                //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                #if DIMENSIONS == 2
                {
                    (*(fluid_beta+0))=flash_vx*cos(ph_phi);
                    (*(fluid_beta+1))=flash_vx*sin(ph_phi);
                    (*(fluid_beta+2))=flash_vy;
                }
                #else
                {
                    (*(fluid_beta+0))=flash_vx;
                    (*(fluid_beta+1))=flash_vy;
                    (*(fluid_beta+2))=flash_vz;
                }
                #endif
                
                /*
                fprintf(fPtr,"FLASH v: %e, %e\n", flash_vx,flash_vy);
                fflush(fPtr);
                */
        
                //fill in photon 4 momentum
                //printf("filling in 4 momentum in photonScatter for photon index %d\n", ph_index);
                if ((ph+ph_index)->type == 's')
                {
                    //printf("The scattering photon is a seed photon w/ comv freq %e Hz.\n", ((ph+ph_index)->comv_p0)*C_LIGHT/PL_CONST);
                    //*nu_c_scatt=((ph+ph_index)->comv_p0)*C_LIGHT/PL_CONST;//dont need this anymore b/c the 's' photon doesnt move from its cell
                
                }
            
            
                *(ph_p+0)=((ph+ph_index)->p0);
                *(ph_p+1)=((ph+ph_index)->p1);
                *(ph_p+2)=((ph+ph_index)->p2);
                *(ph_p+3)=((ph+ph_index)->p3);
                
                *(ph_p_comov+0)=((ph+ph_index)->comv_p0);
                *(ph_p_comov+1)=((ph+ph_index)->comv_p1);
                *(ph_p_comov+2)=((ph+ph_index)->comv_p2);
                *(ph_p_comov+3)=((ph+ph_index)->comv_p3);
            
                //fill in stokes parameters
                *(s+0)=((ph+ph_index)->s0); //I ==1
                *(s+1)=((ph+ph_index)->s1); //Q/I
                *(s+2)=((ph+ph_index)->s2); //U/I
                *(s+3)=((ph+ph_index)->s3); //V/I
        
                /*
                fprintf(fPtr,"Unscattered Photon in Lab frame: %e, %e, %e,%e, %e, %e, %e\nStokes params %e %e %e %e\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3), (ph->r0), (ph->r1), (ph->r2), *(s+0), *(s+1), *(s+2), *(s+3));
                fflush(fPtr);
                fprintf(fPtr,"Fluid Beta: %e, %e, %e\n", *(fluid_beta+0),*(fluid_beta+1), *(fluid_beta+2));
                fflush(fPtr);
                */
        
                //first we bring the photon to the fluid's comoving frame
                //lorentzBoost(fluid_beta, ph_p, ph_p_comov, 'p', fPtr);
                *(ph_p_comov+0)=((ph+ph_index)->comv_p0);
                *(ph_p_comov+1)=((ph+ph_index)->comv_p1);
                *(ph_p_comov+2)=((ph+ph_index)->comv_p2);
                *(ph_p_comov+3)=((ph+ph_index)->comv_p3);
                
                /*
                fprintf(fPtr,"Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
                fflush(fPtr);
         
                fprintf(fPtr, "Before Scattering, In Comov_frame:\n");
                fflush(fPtr);
                fprintf(fPtr, "ph_comov: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                fflush(fPtr);
                */
                        
                //fprintf(fPtr, "Theta: %e Phi %e Lab: x_tilde: %e, %e, %e, y_tilde: %e %e %e\n", theta, phi, *(x_tilde+0), *(x_tilde+1), *(x_tilde+2), *(y_tilde+0), *(y_tilde+1), *(y_tilde+2));
            
                //then rotate the stokes plane by some angle such that we are in the stokes coordinat eystsem after the lorentz boost
                //if (STOKES_SWITCH != 0)
                #if STOKES_SWITCH == ON
                {
                    stokesRotation(fluid_beta, (ph_p+1), (ph_p_comov+1), s, fPtr);
                }
                #endif
                
                //exit(0);
                //second we generate a thermal electron at the correct temperature
                singleElectron(el_p_comov, fluid_temp, ph_p_comov, rand, fPtr);
        
                //fprintf(fPtr,"el_comov: %e, %e, %e,%e\n", *(el_p_comov+0), *(el_p_comov+1), *(el_p_comov+2), *(el_p_comov+3));
                //fflush(fPtr);
         
        
                //third we perform the scattering and save scattered photon 4 monetum in ph_p_comov @ end of function
                event_did_occur=singleScatter(el_p_comov, ph_p_comov, s, rand, fPtr);
            
        
                //fprintf(fPtr,"After Scattering, After Lorentz Boost to Comov frame: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                //fflush(fPtr);
                //event_did_occur=0;
                if (event_did_occur==1)
                {
                    //fprintf(fPtr,"Within the if!\n");
                    //fflush(fPtr);
                
                    //if the scattering occured have to uodate the phtoon 4 momentum. if photon didnt scatter nothing changes
                    //fourth we bring the photon back to the lab frame
                    *(negative_fluid_beta+0)=-1*( *(fluid_beta+0));
                    *(negative_fluid_beta+1)=-1*( *(fluid_beta+1));
                    *(negative_fluid_beta+2)=-1*( *(fluid_beta+2));
                    lorentzBoost(negative_fluid_beta, ph_p_comov, ph_p, 'p',  fPtr);
                    //fprintf(fPtr,"Scattered Photon in Lab frame: %e, %e, %e,%e\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3));
                    //fflush(fPtr);
                    
                    //if (STOKES_SWITCH != 0)
                    #if STOKES_SWITCH == ON
                    {
                        stokesRotation(negative_fluid_beta, (ph_p_comov+1), (ph_p+1), s, fPtr); //rotate to boost back to lab frame
                        
                        //save stokes parameters
                        ((ph+ph_index)->s0)= *(s+0); //I ==1
                        ((ph+ph_index)->s1)= *(s+1);
                        ((ph+ph_index)->s2)= *(s+2);
                        ((ph+ph_index)->s3)= *(s+3);
                    }
                    #endif
                

                    if (((*(ph_p+0))*C_LIGHT/1.6e-9) > 1e4)
                    {
                        fprintf(fPtr,"Extremely High Photon Energy!!!!!!!!\n");
                        fflush(fPtr);
                    }

                    //fprintf(fPtr,"Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
                    //fprintf(fPtr, "Old: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
        
                    //assign the photon its new lab 4 momentum
                    ((ph+ph_index)->p0)=(*(ph_p+0));
                    ((ph+ph_index)->p1)=(*(ph_p+1));
                    ((ph+ph_index)->p2)=(*(ph_p+2));
                    ((ph+ph_index)->p3)=(*(ph_p+3));
                    
                    //assign it the comoving frame 4 momentum
                    ((ph+ph_index)->comv_p0)=(*(ph_p_comov+0));
                    ((ph+ph_index)->comv_p1)=(*(ph_p_comov+1));
                    ((ph+ph_index)->comv_p2)=(*(ph_p_comov+2));
                    ((ph+ph_index)->comv_p3)=(*(ph_p_comov+3));
                    
                    //printf("Done assigning values to original struct\n");
        
                    //incremement that photons number of scatterings
                    ((ph+ph_index)->num_scatt)+=1;
                    *frame_scatt_cnt+=1; //incrememnt total number of scatterings
                    
                    
                    
                
                }
                /*
            }
            else
            {
                //the photon will be absorbed
                event_did_occur=1;
                *frame_abs_cnt+=1;
                //set weight=0 (to make sure it cant affect anything if it somehow gets saved) and for us to knwo which photons have been absorbed
                ((ph+ph_index)->weight)=0;
                ((ph+ph_index)->nearest_block_index)=-1;
                fprintf(fPtr, "Photon %d In the absorption part of if-else.\n", ph_index);
                //exit(0);
            } wont need this since we are not dealing with individual absorption events
                 */
        }
        else
        {
            // if the photon scatt_time > dt_max
            //have to adjust the time properly so that the time si now appropriate for the next frame
            scatt_time=dt_max;
            updatePhotonPosition(ph, num_ph, scatt_time-old_scatt_time, fPtr); 
            event_did_occur=1; //set equal to 1 to get out of the loop b/c other subsequent photons will have scatt_time > dt_max
            
        }
    
        old_scatt_time=scatt_time;
        i++;
	}
    //exit(0);
    *scattered_ph_index=ph_index; //save the index of the photon that was scattered
    
    //fprintf(fPtr,"scattered_ph_index: %d %d\n", *scattered_ph_index, (*(sorted_indexes+i-1)));
    //fflush(fPtr);
    
    free(el_p_comov); 
    free(ph_p_comov);
    free(fluid_beta); 
    free(negative_fluid_beta);
    free(ph_p);
    free(s);
    ph_p=NULL;negative_fluid_beta=NULL;ph_p_comov=NULL; el_p_comov=NULL;
    
    //retrun total time elapsed to scatter a photon
    return scatt_time;
}

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr)
{
    //generates an electron with random energy 
    double factor=0, gamma=0;
    double y_dum=0, f_x_dum=0, x_dum=0, beta_x_dum=0, beta=0, phi=0, theta=0, ph_theta=0, ph_phi=0;
    gsl_matrix *rot= gsl_matrix_calloc (3, 3); //create matrix thats 3x3 to do rotation 
    gsl_vector_view el_p_prime ; //create vector to hold rotated electron 4 momentum
    gsl_vector *result=gsl_vector_alloc (3);
    
    //fprintf(fPtr, "Temp in singleElectron: %e\n", temp);
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
            
            f_x_dum=pow(x_dum,2)*(beta_x_dum/gsl_sf_bessel_Kn (2, 1.0/factor))*exp(-1*x_dum/factor); //
            //fprintf(fPtr,"Choosing a Gamma: xdum: %e, f_x_dum: %e, y_dum: %e\n", x_dum, f_x_dum, y_dum);
        }
        gamma=x_dum;
        
    }
    else
    {

        //printf("In else\n");
        factor=pow(K_B*temp/M_EL,0.5);
        //calculate a random gamma from 3 random velocities drawn from a gaussian distribution with std deviation of "factor"
        gamma=pow( 1- (pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+ pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)  ) ,-0.5); //each vel direction is normal distribution -> maxwellian when multiplied
    }
    
    //fprintf(fPtr,"Chosen Gamma: %e\n",gamma);
    
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
    //fprintf(fPtr,"Beta: %e\tPhi: %e\tTheta: %e\n",beta,phi, theta);
    //fill in electron 4 momentum NOT SURE WHY THE ORDER IS AS SUCH SEEMS TO BE E/c, pz,py,px!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    *(el_p+0)=gamma*(M_EL)*(C_LIGHT);
    *(el_p+1)=gamma*(M_EL)*(C_LIGHT)*beta*cos(theta);
    *(el_p+2)=gamma*(M_EL)*(C_LIGHT)*beta*sin(theta)*sin(phi);
    *(el_p+3)=gamma*(M_EL)*(C_LIGHT)*beta*sin(theta)*cos(phi);
    
    //printf("Old: %e, %e, %e,%e\n", *(el_p+0), *(el_p+1), *(el_p+2), *(el_p+3));
    
    el_p_prime=gsl_vector_view_array((el_p+1), 3);
    
    //find angles of photon NOT SURE WHY WERE CHANGING REFERENCE FRAMES HERE???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ph_phi=atan2(*(ph_p+2), *(ph_p+3)); //Double Check
    ph_theta=atan2(pow( pow(*(ph_p+2),2)+  pow(*(ph_p+3),2) , 0.5) , (*(ph_p+1)) );
    
    //printf("Calculated Photon phi and theta in singleElectron:%e, %e\n", ph_phi, ph_theta);
    
    //fill in rotation matrix to rotate around x axis to get rid of phi angle
    gsl_matrix_set(rot, 1,1,1);
    gsl_matrix_set(rot, 2,2,cos(ph_theta));
    gsl_matrix_set(rot, 0,0,cos(ph_theta));
    gsl_matrix_set(rot, 0,2,-sin(ph_theta));
    gsl_matrix_set(rot, 2,0,sin(ph_theta));
    gsl_blas_dgemv(CblasNoTrans, 1, rot, &el_p_prime.vector, 0, result);
    
    /*
    printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot, 0,0), gsl_matrix_get(rot, 0,1), gsl_matrix_get(rot, 0,2));
    printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot, 1,0), gsl_matrix_get(rot, 1,1), gsl_matrix_get(rot, 1,2));
    printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot, 2,0), gsl_matrix_get(rot, 2,1), gsl_matrix_get(rot, 2,2));

    printf("Middle: %e, %e, %e,%e\n", *(el_p+0), gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2));
    */
    
    gsl_matrix_set_all(rot,0);
    
    gsl_matrix_set(rot, 0,0,1);
    gsl_matrix_set(rot, 1,1,cos(-ph_phi));
    gsl_matrix_set(rot, 2,2,cos(-ph_phi));
    gsl_matrix_set(rot, 1,2,-sin(-ph_phi));
    gsl_matrix_set(rot, 2,1,sin(-ph_phi));
    gsl_blas_dgemv(CblasNoTrans, 1, rot, result, 0, &el_p_prime.vector);
    /*
    printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot, 0,0), gsl_matrix_get(rot, 0,1), gsl_matrix_get(rot, 0,2));
    printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot, 1,0), gsl_matrix_get(rot, 1,1), gsl_matrix_get(rot, 1,2));
    printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot, 2,0), gsl_matrix_get(rot, 2,1), gsl_matrix_get(rot, 2,2));
    printf("Final EL_P_vec: %e, %e, %e,%e\n", *(el_p+0), gsl_vector_get(&el_p_prime.vector,0), gsl_vector_get(&el_p_prime.vector,1), gsl_vector_get(&el_p_prime.vector,2));
    */
    
    
    gsl_matrix_free (rot);gsl_vector_free(result);
}


int singleScatter(double *el_comov, double *ph_comov, double *s, gsl_rng * rand, FILE *fPtr)
{
    //This routine performs a scattering between a photon and a moving electron.
    int i=0, scattering_occured=0;
    double dotprod_1; //to test orthogonality
    double *z_axis_electron_rest_frame=malloc(3*sizeof(double)); //basis vector of the z axis in the elctron rest frame
    double *el_v=malloc(3*sizeof(double));
    double *negative_el_v=malloc(3*sizeof(double));
    double *ph_p_prime=malloc(4*sizeof(double));//use this to keep track of how the ph 4 momentum changes with each rotation
    double *el_p_prime=malloc(4*sizeof(double));
    double phi0=0, phi1=0, phi=0, theta=0;
    double y_dum, f_x_dum, x_dum;
    double x_tilde[3]={0,0,0}, y_tilde[3]={0,0,0}, x_tilde_new[3]={0,0,0}, y_tilde_new[3]={0,0,0};//initalize arrays to hold stokes coordinate system
    gsl_matrix *rot0= gsl_matrix_calloc (3, 3); //create matricies thats 3x3 to do rotations
    gsl_matrix *rot1= gsl_matrix_calloc (3, 3);
    gsl_matrix *scatt= gsl_matrix_calloc (4, 4); //fano's matrix for scattering stokes parameters
    gsl_vector *scatt_result=gsl_vector_alloc (4);
    gsl_vector *result0=gsl_vector_alloc (3); //vectors to hold results of rotations
    gsl_vector *result1=gsl_vector_alloc (3);
    gsl_vector *result=gsl_vector_alloc (4);
    gsl_vector *whole_ph_p=gsl_vector_alloc (4);
    gsl_vector *ph_p_orig=gsl_vector_alloc (4) ;//vector to hold the original incoming photon velocity vector in the electron rest frame
    gsl_vector_view ph_p ;//create vector to hold comoving photon and electron 4 momentum
    gsl_vector_view el_p ;
    gsl_vector_view stokes, test, test_x, test_y;
    /*
     Dont need these vectors anymore, plus didnt have code to free allocations so it was causing memory leaks
    gsl_vector *result0_x=gsl_vector_alloc (3); //vectors to hold results of rotations for stokes coordinates
    gsl_vector *result1_x=gsl_vector_alloc (3);
    gsl_vector *result0_y=gsl_vector_alloc (3); //vectors to hold results of rotations for stokes coordinates
    gsl_vector *result1_y=gsl_vector_alloc (3);
     */
    
    //fill in z-axis basis vector
    *(z_axis_electron_rest_frame+0)=0;
    *(z_axis_electron_rest_frame+1)=0;
    *(z_axis_electron_rest_frame+2)=1;
    
    /* was for testing against Kraw
    *(s+0)=1; //should be 1.0
    *(s+1)=1;
    *(s+2)=0;
    *(s+3)=0;
    
    *(ph_comov+0)=PL_CONST*1e12/C_LIGHT;
    *(ph_comov+1)=0; //set values of photon prime momentum from doing the scattering to use the vector view of it in dot product
    *(ph_comov+2)=0;
    *(ph_comov+3)=PL_CONST*1e12/C_LIGHT;
    
    theta=85*M_PI/180;
    phi=0;
    dotprod_1=pow(1-(pow(100, -2.0)) ,0.5);
    *(el_comov+0)=100*M_EL*C_LIGHT;
    *(el_comov+1)=100*M_EL*C_LIGHT*dotprod_1*sin(theta)*cos(phi); //set values of photon prime momentum from doing the scattering to use the vector view of it in dot product
    *(el_comov+2)=100*M_EL*C_LIGHT*dotprod_1*sin(theta)*sin(phi);
    *(el_comov+3)=100*M_EL*C_LIGHT*dotprod_1*cos(theta);
     */
    
    //fill in electron velocity array and photon 4 momentum
    *(el_v+0)=(*(el_comov+1))/(*(el_comov+0));
    *(el_v+1)=(*(el_comov+2))/(*(el_comov+0));
    *(el_v+2)=(*(el_comov+3))/(*(el_comov+0));
    //printf("el_v: %e, %e, %e\n", *(el_v+0), *(el_v+1), *(el_v+2));
    
    //lorentz boost into frame where the electron is stationary
    lorentzBoost(el_v, el_comov, el_p_prime, 'e', fPtr);
    lorentzBoost(el_v, ph_comov, ph_p_prime, 'p', fPtr);
    //printf("New ph_p in electron rest frame: %e, %e, %e,%e\n", *(ph_p_prime+0), *(ph_p_prime+1), *(ph_p_prime+2), *(ph_p_prime+3));
    
    //rotate 'stokes plane'
    //if (STOKES_SWITCH != 0)
    #if STOKES_SWITCH == ON
    {
        i=1;
        stokesRotation(el_v, (ph_comov+1), (ph_p_prime+1), s, fPtr);
        stokes=gsl_vector_view_array(s, 4);
        i=0;
    }
    #endif
    
    //printf(fPtr, "y_tilde: %e, %e, %e\n", *(y_tilde+0), *(y_tilde+1), *(y_tilde+2));
    
    ph_p=gsl_vector_view_array((ph_p_prime+1), 3);
    el_p=gsl_vector_view_array(el_p_prime,4);
    
    gsl_vector_set(ph_p_orig, 0, *(ph_p_prime+0));
    gsl_vector_set(ph_p_orig, 1, *(ph_p_prime+1));
    gsl_vector_set(ph_p_orig, 2, *(ph_p_prime+2));
    gsl_vector_set(ph_p_orig, 3, *(ph_p_prime+3));
    
    //gsl_blas_ddot(&y_tilde_rot.vector, &ph_p.vector, &dotprod_1);
    //fprintf(fPtr, "After lorentz boost Angle between the  y_tilde_rot and the photon velocity vector is: %e\n", acos(dotprod_1/ gsl_blas_dnrm2(&ph_p.vector))*180/M_PI);
    
    phi0=atan2(*(ph_p_prime+2), *(ph_p_prime+1) );
    //fprintf(fPtr,"Photon Phi: %e\n", phi0);
    //rotate the axes so that the photon incomes along the x-axis
    gsl_matrix_set(rot0, 2,2,1);
    gsl_matrix_set(rot0, 0,0,cos(-phi0));
    gsl_matrix_set(rot0, 1,1,cos(-phi0));
    gsl_matrix_set(rot0, 0,1,-sin(-phi0));
    gsl_matrix_set(rot0, 1,0,sin(-phi0));
    gsl_blas_dgemv(CblasNoTrans, 1, rot0, &ph_p.vector, 0, result0);
    
    //printf("Before Scatter rot0: stokes x=(%e, %e, %e) y=(%e, %e, %e)", gsl_vector_get(result0_x,0), gsl_vector_get(result0_x,1), gsl_vector_get(result0_x,2), gsl_vector_get(result0_y,0), gsl_vector_get(result0_y,1), gsl_vector_get(result0_y,2));
    
    //fprintf(fPtr, "y_tilde: %e, %e, %e y_tilde_rot_result: %e, %e, %e\n", *(y_tilde+0), *(y_tilde+1), *(y_tilde+2), gsl_vector_get(y_tilde_rot_result,0), gsl_vector_get(y_tilde_rot_result,1), gsl_vector_get(y_tilde_rot_result,2));
    
    /*
     printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot0, 0,0), gsl_matrix_get(rot0, 0,1), gsl_matrix_get(rot0, 0,2));
     printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot0, 1,0), gsl_matrix_get(rot0, 1,1), gsl_matrix_get(rot0, 1,2));
     printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot0, 2,0), gsl_matrix_get(rot0, 2,1), gsl_matrix_get(rot0, 2,2));
     */
    
    //set values of ph_p_prime equal to the result and get new phi from result
    *(ph_p_prime+1)=gsl_vector_get(result0,0);
    *(ph_p_prime+2)=0;//gsl_vector_get(result,1); //just directly setting it to 0 now?
    *(ph_p_prime+3)=gsl_vector_get(result0,2);
    
    phi1=atan2(gsl_vector_get(result0,2), gsl_vector_get(result0,0));
    
    
    //printf("rotation 1: %e, %e, %e\n",  *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    //fprintf(fPtr, "Photon Phi: %e\n", phi1);
    //printf("make sure the vector view is good: %e, %e, %e,%e\n", *(ph_p_prime+0), gsl_vector_get(&ph_p.vector,0), gsl_vector_get(&ph_p.vector,1), gsl_vector_get(&ph_p.vector,2));
    
    
    //rotate around y to bring it all along x
    gsl_matrix_set(rot1, 1,1,1);
    gsl_matrix_set(rot1, 0,0,cos(-phi1));
    gsl_matrix_set(rot1, 2,2,cos(-phi1));
    gsl_matrix_set(rot1, 0,2,-sin(-phi1));
    gsl_matrix_set(rot1, 2,0,sin(-phi1));
    gsl_blas_dgemv(CblasNoTrans, 1, rot1, &ph_p.vector, 0, result1);
    
    //fprintf(fPtr, "y_tilde: %e, %e, %e y_tilde_rot vector view: %e, %e, %e\n", *(y_tilde+0), *(y_tilde+1), *(y_tilde+2), gsl_vector_get(&y_tilde_rot.vector,0), gsl_vector_get(&y_tilde_rot.vector,1), gsl_vector_get(&y_tilde_rot.vector,2));
    
    /*
     printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot1, 0,0), gsl_matrix_get(rot1, 0,1), gsl_matrix_get(rot1, 0,2));
     printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot1, 1,0), gsl_matrix_get(rot1, 1,1), gsl_matrix_get(rot1, 1,2));
     printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot1, 2,0), gsl_matrix_get(rot1, 2,1), gsl_matrix_get(rot1, 2,2));
     */
    
    //set values of ph_p_prime equal to the result and get new phi from result
    *(ph_p_prime+1)=*(ph_p_prime+0);//why setting it to the energy?
    *(ph_p_prime+2)=gsl_vector_get(result1,1);
    *(ph_p_prime+3)=0; //just directly setting it to 0 now?
    
    //printf("rotation 2: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    
    //know that the stokes y axis is in -y_hat direction and stokes x asis is in the z_hat direction due to rotations and making inclimg photn come along x_hat direction, dont need to rotate the stokes plane/vector. this happens as the rotations occur (tested in python code)
    //double checking here
    //printf("Before Scatter: stokes x=(%e, %e, %e) y=(%e, %e, %e) ph_p=(%e, %e, %e, %e)\n", gsl_vector_get(result1_x,0), gsl_vector_get(result1_x,1), gsl_vector_get(result1_x,2), gsl_vector_get(result1_y,0), gsl_vector_get(result1_y,1), gsl_vector_get(result1_y,2), *(ph_p_prime+0), *(ph_p_prime+1), *(ph_p_prime+2), *(ph_p_prime+3));
    
    
    //determine if the scattering will occur between photon and electron
    //scattering_occured=comptonScatter(&theta, &phi, rand, fPtr); //determine the angles phi and theta for the photon to scatter into using thompson differential cross section
    scattering_occured=kleinNishinaScatter(&theta, &phi, *(ph_p_prime+0), *(s+1), *(s+2), rand, fPtr);//determine the angles phi and theta for the photon to scatter into using KN differential cross section, if the photon will end up scattering
    
    //fprintf(fPtr,"Phi: %e, Theta: %e\n", phi, theta);
    //theta=2.4475668271885342;
    //phi=4.014719957630734;
    //*(s+0)=1; //should be 1.0
    //*(s+1)=1;
    //*(s+2)=0;
    //*(s+3)=0;
    
    
    if (scattering_occured==1)
    {
        //perform scattering and compute new 4-momenta of electron and photon
        //scattered photon 4 momentum
        gsl_vector_set(result, 0, (*(ph_p_prime+0))/(1+ (( (*(ph_p_prime+0))*(1-cos(theta)) )/(M_EL*C_LIGHT )) ) ); // scattered energy of photon
        gsl_vector_set(result, 1, gsl_vector_get(result,0)*cos(theta) );
        gsl_vector_set(result, 2, gsl_vector_get(result,0)*sin(theta)*sin(phi) );//assume phi is clockwise from z to y
        gsl_vector_set(result, 3, gsl_vector_get(result,0)*sin(theta)*cos(phi) );
        //fprintf(fPtr, "New ph_p0=%e Old= %e\n", gsl_vector_get(result,0), *(ph_p_prime+0));
        //gsl_vector_fprintf(fPtr,result, "%e" );
        
        //recalc x_tilde from rotation about y by angle theta do x_tilde=y_tilde X v_ph
        //test =gsl_vector_view_array(gsl_vector_ptr(result, 1), 3);
        
        //scatt_result is a dummy, dont need to change the stokes parameters here, just need to find the axis such that y is out of the plane of k_o-k see Ito figure 12 in polarized emission from stratisfied jets
        
        //gsl_blas_ddot(&y_tilde_rot.vector, &test.vector, &dotprod_1);
        //fprintf(fPtr, "Angle between the  y_tilde_rot and the photon velocity vector is: %e\n", acos(dotprod_1/ gsl_blas_dnrm2(&test.vector))*180/M_PI);
        //gsl_vector_fprintf(fPtr,&y_tilde_rot.vector, "%e" );
        //gsl_vector_fprintf(fPtr,&x_tilde_rot.vector, "%e" );
        
        //exit(0);
        //calculate electron 4 momentum
        //prescattered photon 4 momentum
        gsl_vector_set(whole_ph_p, 0, (*(ph_p_prime+0)));
        gsl_vector_set(whole_ph_p, 1, (*(ph_p_prime+1)));
        gsl_vector_set(whole_ph_p, 2, (*(ph_p_prime+2)));
        gsl_vector_set(whole_ph_p, 3, (*(ph_p_prime+3)));
        
        gsl_vector_sub(whole_ph_p,result); //resut is saved into ph_p vector, unscattered-scattered 4 mometum of photon
        gsl_vector_add(&el_p.vector ,whole_ph_p);
        /*
         printf("After scattering:\n");
         printf("el_p: %e, %e, %e,%e\n", gsl_vector_get(&el_p.vector,0), gsl_vector_get(&el_p.vector,1), gsl_vector_get(&el_p.vector,2), gsl_vector_get(&el_p.vector,3));
         printf("ph_p: %e, %e, %e,%e\n", gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2), gsl_vector_get(result,3));
         */
        
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
        /*
         printf("Photon Phi: %e\n", phi1);
         printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot1, 0,0), gsl_matrix_get(rot1, 0,1), gsl_matrix_get(rot1, 0,2));
         printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot1, 1,0), gsl_matrix_get(rot1, 1,1), gsl_matrix_get(rot1, 1,2));
         printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot1, 2,0), gsl_matrix_get(rot1, 2,1), gsl_matrix_get(rot1, 2,2));
         */
        
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
        
        /*
         printf("Photon Phi: %e\n", phi0);
         printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot0, 0,0), gsl_matrix_get(rot0, 0,1), gsl_matrix_get(rot0, 0,2));
         printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot0, 1,0), gsl_matrix_get(rot0, 1,1), gsl_matrix_get(rot0, 1,2));
         printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot0, 2,0), gsl_matrix_get(rot0, 2,1), gsl_matrix_get(rot0, 2,2));
         */
        
        //do the scattering of the stokes vector
        //rotate it by phi and then scatter it and rotate back and then renormalize it such that i=1
        //if (STOKES_SWITCH != 0)
        #if STOKES_SWITCH == ON
        {
            //orient the stokes coordinate system such that its perpendicular to the scattering plane
            findXY(gsl_vector_ptr(ph_p_orig, 1),z_axis_electron_rest_frame, x_tilde, y_tilde);
            findXY(gsl_vector_ptr(result0,0),gsl_vector_ptr(ph_p_orig, 1), x_tilde_new, y_tilde_new);
            phi=findPhi(x_tilde, y_tilde, x_tilde_new, y_tilde_new);
            mullerMatrixRotation(phi, s, fPtr);
            
            //find the theta between the incoming and scattered photons, by doing dot product and taking arccos of it
            theta=acos((gsl_vector_get(ph_p_orig,1)*gsl_vector_get(result0,0)+gsl_vector_get(ph_p_orig,2)*gsl_vector_get(result0,1)+gsl_vector_get(ph_p_orig,3)*gsl_vector_get(result0,2) )/(gsl_vector_get(ph_p_orig,0)*(*(ph_p_prime+0)))  );
            
            //do the scattering of the stokes parameters
            gsl_matrix_set(scatt, 0,0,1.0+pow(cos(theta), 2.0)+((1-cos(theta))*(gsl_vector_get(ph_p_orig,0) - gsl_vector_get(result,0))/(M_EL*C_LIGHT ) ) ); //following lundman's matrix
            gsl_matrix_set(scatt, 0,1, sin(theta)*sin(theta));
            gsl_matrix_set(scatt, 1,0, sin(theta)*sin(theta));
            gsl_matrix_set(scatt, 1,1,1.0+cos(theta)*cos(theta));
            gsl_matrix_set(scatt, 2,2, 2.0*cos(theta));
            gsl_matrix_set(scatt, 3,3, 2.0*cos(theta)+ ((cos(theta))*(1-cos(theta))*(gsl_vector_get(ph_p_orig,0) - gsl_vector_get(result,0))/(M_EL*C_LIGHT )) );
            //gsl_matrix_scale(scatt, (gsl_vector_get(result,0)/(*(ph_p_prime+0)))*((gsl_vector_get(result,0)/(*(ph_p_prime+0))))*0.5*3*THOM_X_SECT/(8*M_PI) ); //scale the matrix by 0.5*r_0^2 (\epsilon/\epsilon_0)^2 DONT NEED THIS BECAUSE WE NORMALIZE STOKES VECTOR SO THIS CANCELS ITSELF OUT
            gsl_blas_dgemv(CblasNoTrans, 1, scatt, &stokes.vector, 0, scatt_result);
            /*
             fprintf(fPtr,"before s: %e, %e, %e,%e\n", gsl_vector_get(&stokes.vector,0), gsl_vector_get(&stokes.vector,1), gsl_vector_get(&stokes.vector,2), gsl_vector_get(&stokes.vector,3));
             fprintf(fPtr,"Scatt Matrix 0: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 0,0), gsl_matrix_get(scatt, 0,1), gsl_matrix_get(scatt, 0,2), gsl_matrix_get(scatt, 0,3));
             fprintf(fPtr,"Scatt Matrix 1: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 1,0), gsl_matrix_get(scatt, 1,1), gsl_matrix_get(scatt, 1,2), gsl_matrix_get(scatt, 1,3));
             fprintf(fPtr,"Scatt Matrix 2: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 2,0), gsl_matrix_get(scatt, 2,1), gsl_matrix_get(scatt, 2,2), gsl_matrix_get(scatt, 2,3));
             fprintf(fPtr,"Scatt Matrix 3: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 3,0), gsl_matrix_get(scatt, 3,1), gsl_matrix_get(scatt, 3,2), gsl_matrix_get(scatt, 3,3));
             fprintf(fPtr,"s: %e, %e, %e,%e\n", gsl_vector_get(scatt_result,0), gsl_vector_get(scatt_result,1), gsl_vector_get(scatt_result,2), gsl_vector_get(scatt_result,3));
             */
            //normalize and rotate back
            *(s+0)=gsl_vector_get(scatt_result,0)/gsl_vector_get(scatt_result,0); //should be 1.0
            *(s+1)=gsl_vector_get(scatt_result,1)/gsl_vector_get(scatt_result,0);
            *(s+2)=gsl_vector_get(scatt_result,2)/gsl_vector_get(scatt_result,0);
            *(s+3)=gsl_vector_get(scatt_result,3)/gsl_vector_get(scatt_result,0);
            //fprintf(fPtr,"s after norm: %e, %e, %e,%e\n", gsl_vector_get(&stokes.vector,0), gsl_vector_get(&stokes.vector,1), gsl_vector_get(&stokes.vector,2), gsl_vector_get(&stokes.vector,3));
            
            //need to find current stokes coordinate system defined in the plane of k-k_0
            findXY(gsl_vector_ptr(result0,0),gsl_vector_ptr(ph_p_orig, 1), x_tilde, y_tilde);
            
            //then find the new coordinate system between scattered photon 4 onetum and the z axis
            findXY(gsl_vector_ptr(result0,0),z_axis_electron_rest_frame, x_tilde_new, y_tilde_new);
            
            //find phi to transform between the two coodinate systems
            phi=findPhi(x_tilde, y_tilde, x_tilde_new, y_tilde_new);
            
            //do the rotation
            mullerMatrixRotation(phi, s, fPtr);
        }
        #endif
        
        //now update the array with the new scattered photon 4 monetum
        *(ph_p_prime+1)=gsl_vector_get(result0,0);
        *(ph_p_prime+2)=gsl_vector_get(result0,1);
        *(ph_p_prime+3)=gsl_vector_get(result0,2);
        
        //gsl_blas_ddot(&y_tilde_rot.vector, &ph_p.vector, &dotprod_1);
        //fprintf(fPtr, "Angle between the  y_tilde_rot and the photon velocity vector is: %e\n", acos(dotprod_1/ gsl_blas_dnrm2(&ph_p.vector))*180/M_PI);
        
        //printf("Undo rotation 1: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
        //deboost photon to lab frame
        *(negative_el_v+0)=(-1*(*(el_v+0)));
        *(negative_el_v+1)=(-1*(*(el_v+1)));
        *(negative_el_v+2)=(-1*(*(el_v+2)));
        
        lorentzBoost(negative_el_v, ph_p_prime, ph_comov, 'p', fPtr);
        //printf("Undo boost 1: %e, %e, %e, %e\n",  *(ph_comov+0), *(ph_comov+1),  *(ph_comov+2),  *(ph_comov+3));
        
        //dont need to find stokes vector and do previosu rotations, can just find the stokes coordinates in function because the stokes coordinate vectors rotate with the photon vector and no rotations to a new stokes coordinate system are needed
        //if (STOKES_SWITCH != 0)
        #if STOKES_SWITCH == ON
        {
            stokesRotation(negative_el_v, (ph_p_prime+1), (ph_comov+1), s, fPtr);
        }
        #endif
        
        //exit(0);
    }
    
    gsl_matrix_free(rot0); gsl_matrix_free(rot1);gsl_matrix_free(scatt);gsl_vector_free(result0);gsl_vector_free(result1);gsl_vector_free(result);
    gsl_vector_free(scatt_result);gsl_vector_free(ph_p_orig);
    gsl_vector_free(whole_ph_p);free(ph_p_prime);free(el_p_prime);free(el_v); free(negative_el_v); free(z_axis_electron_rest_frame);
    
    return scattering_occured;
}

int comptonScatter(double *theta, double *phi, gsl_rng * rand, FILE *fPtr)
{
    
        double y_dum, f_x_dum, x_dum;
        
        //generate random theta and phi angles for scattering
        *phi=gsl_rng_uniform(rand)*2*M_PI;
        //printf("Phi: %e\n", phi);
    
        y_dum=1; //initalize loop to get a random theta
        f_x_dum=0;
        while (y_dum>f_x_dum)
        {
            y_dum=gsl_rng_uniform(rand)*1.09;
            x_dum=gsl_rng_uniform(rand)*M_PI;
            f_x_dum=sin(x_dum)*(1+pow(cos(x_dum),2));
        }
        *theta=x_dum;
        
        return 1;
}


int kleinNishinaScatter(double *theta, double *phi, double p0, double q, double u, gsl_rng * rand, FILE *fPtr)
{
    //sample theta using:  https://doi.org/10.13182/NSE11-57
    double phi_dum=0, cos_theta_dum=0, f_phi_dum=0, f_cos_theta_dum=0, f_theta_dum=0, phi_y_dum=0, cos_theta_y_dum=0, KN_x_section_over_thomson_x_section=0, rand_num=0;
    double mu=0, phi_norm=0, phi_max=0, norm=0;
    int will_scatter=0;
    double energy_ratio=  p0/(M_EL*C_LIGHT ); //h*nu / mc^2 , units of p0 is erg/c 
    
    //determine the KN cross section over the thomson cross section From RYBICKI AND LIGHTMAN pg 197
    KN_x_section_over_thomson_x_section= (3.0/4.0)*(  (  ((1+energy_ratio)/ pow(energy_ratio,3.0))*(((2*energy_ratio)*(1+energy_ratio)/(1+2*energy_ratio)) - log(1+2*energy_ratio)))  + (log(1+2*energy_ratio)/(2*energy_ratio)) - ((1+3*energy_ratio)/pow((1+2*energy_ratio),2.0))  );
    rand_num=gsl_rng_uniform(rand);
    
    //fprintf(fPtr,"Rand: %e, p0: %e, X: %e, Ratio: %e\n", rand_num, p0*C_LIGHT, energy_ratio, KN_x_section_over_thomson_x_section);
    //fflush(fPtr);
    
    if ((rand_num<= KN_x_section_over_thomson_x_section) || (p0 < 1e-2*(M_EL*C_LIGHT ) ))
    {
        //include last condition so low energy seed phtoons can scatter (as they should under thompson scattering), calculating KN_x_section_over_thomson_x_section incurs numerical error at very low frequencies
        //fprintf(fPtr,"In If!\n");
        //fflush(fPtr);
    
        //sample a theta and phi from the differential cross sections
        phi_y_dum=1; //initalize loop to get a random phi and theta
        cos_theta_y_dum=1;
        f_cos_theta_dum=0;
        f_phi_dum=0;
        
        while ((cos_theta_y_dum>f_cos_theta_dum))
        {
            //do phi and theta seperately, sample theta using:  https://doi.org/10.13182/NSE11-57
            cos_theta_y_dum=gsl_rng_uniform(rand)*2;
            cos_theta_dum=gsl_rng_uniform(rand)*2-1;
            f_cos_theta_dum=pow((1+energy_ratio*(1-cos_theta_dum)),-2)*(energy_ratio*(1-cos_theta_dum)+(1/(1+energy_ratio*(1-cos_theta_dum))) + cos_theta_dum*cos_theta_dum);

            //fprintf(fPtr,"theta_y_dum: %e, theta_dum: %e, mu: %e, f_theta_dum: %e\n", theta_y_dum, theta_dum, mu, f_theta_dum);
            //fflush(fPtr);
        }
        *theta=acos(cos_theta_dum);
        mu=1+energy_ratio*(1-cos(*theta));
        f_theta_dum=(pow(mu, -1.0) + pow(mu, -3.0) - pow(mu, -2.0)*pow(sin(*theta), 2.0))*sin(*theta);
        
        while ((phi_y_dum>f_phi_dum) )
        {
            
            #if STOKES_SWITCH == OFF
            {
                //not considering polarization therefore can jjst sample between 0 and 2*pi evenly
                phi_dum=gsl_rng_uniform(rand)*2*M_PI;
                phi_y_dum=-1; // this is to exit the while statement
                
                //fprintf(fPtr," phi_dum: %e\n", phi_dum);
                //fflush(fPtr);

            }
            #else
            {
                if (u!=0 && q!=0)
                {
                    //if we are considering polarization calulate the norm for the distributiion to be between 1 and 0
                    phi_max=abs(atan2(-u,q))/2.0;
                    norm=(f_theta_dum + pow(mu, -2.0)*pow(sin(*theta), 3.0) * (q*cos(2*phi_max)-u*sin(2*phi_max)));
                    //fprintf(fPtr,"norm: %e\n", norm);
                    //fflush(fPtr);
                    
                    phi_y_dum=gsl_rng_uniform(rand);
                    phi_dum=gsl_rng_uniform(rand)*2*M_PI;
                    f_phi_dum=(f_theta_dum + pow(mu, -2.0)*pow(sin(*theta), 3.0) * (q*cos(2*phi_dum)-u*sin(2*phi_dum)))/norm; //signs on q and u based on Lundman/ McMaster
                    
                    //fprintf(fPtr,"phi_y_dum: %e, theta_dum: %e, mu: %e, f_theta_dum: %e, phi_dum: %e, f_phi_dum: %e, u: %e, q: %e\n", phi_y_dum, theta_dum, mu, f_theta_dum, phi_dum, f_phi_dum, u, q);
                    //fflush(fPtr);
                }
                else
                {
                    phi_dum=gsl_rng_uniform(rand)*2*M_PI;
                    phi_y_dum=-1; // this is to exit the while statement

                }
            }
            #endif
            
        }
        *phi=phi_dum;
        
        will_scatter=1;
    }
    else
    {
        will_scatter=0;
    }
    
    return will_scatter;
}

double averagePhotonEnergy(struct photon *ph, int num_ph)
{
    //to calculate weighted photon energy in ergs
    int i=0, num_thread=omp_get_num_threads();
    double e_sum=0, w_sum=0;
    
    #pragma omp parallel for reduction(+:e_sum) reduction(+:w_sum)
    for (i=0;i<num_ph;i++)
    {
        if (((ph+i)->weight != 0) && ((ph+i)->nearest_block_index != -1))
        {
            e_sum+=(((ph+i)->p0)*((ph+i)->weight));
            w_sum+=((ph+i)->weight);
        }
    }
    
    return (e_sum*C_LIGHT)/w_sum;
}

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg, double *r_avg  )
{
    int temp_max=0, temp_min=INT_MAX,  i=0, num_thread=omp_get_num_threads();
    double sum=0, avg_r_sum=0;
    
    //printf("Num threads: %d", num_thread);
    #pragma omp parallel for num_threads(num_thread) reduction(min:temp_min) reduction(max:temp_max) reduction(+:sum) reduction(+:avg_r_sum)
    for (i=0;i<ph_num;i++)
    {
        if (((ph+i)->weight != 0) && ((ph+i)->nearest_block_index != -1))
        {
            sum+=((ph+i)->num_scatt);
            avg_r_sum+=pow(((ph+i)->r0)*((ph+i)->r0) + ((ph+i)->r1)*((ph+i)->r1) + ((ph+i)->r2)*((ph+i)->r2), 0.5);
            
            if (((ph+i)->num_scatt) > temp_max )
            {
                temp_max=((ph+i)->num_scatt);
                //printf("The new max is: %d\n", temp_max);
            }
            
            //if ((i==0) || (((ph+i)->num_scatt)<temp_min))
            if (((ph+i)->num_scatt)<temp_min)
            {
                temp_min=((ph+i)->num_scatt);
                //printf("The new min is: %d\n", temp_min);
            }
        }
        
    }
    //printf("The  min outside the loop is: %d\n", temp_min);
    //exit(0);
    
    *avg=sum/ph_num;
    *r_avg=avg_r_sum/ph_num;
    *max=temp_max;
    *min=temp_min;
    
}

void cylindricalPrep(double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array)
{
    double  gamma_infinity=100, t_comov=1*pow(10, 5), ddensity=3e-7;// the comoving temperature in Kelvin, and the comoving density in g/cm^2
    int i=0;
    double vel=pow(1-pow(gamma_infinity, -2.0) ,0.5), lab_dens=gamma_infinity*ddensity;
    
    for (i=0; i<num_array;i++)
    {
        *(gamma+i)=gamma_infinity;
        *(vx+i)=0;
        *(vy+i)=vel;
        *(dens+i)=ddensity;
        *(dens_lab+i)=lab_dens;
        *(pres+i)=(A_RAD*pow(t_comov, 4.0))/(3*pow(C_LIGHT, 2.0)); 
        *(temp+i)=pow(3*(*(pres+i))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0); //just assign t_comov
    }
    
}

void sphericalPrep(double *r,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr)
{
    double  gamma_infinity=100, lumi=1e56, r00=1e8; //shopuld be 10^57
    double vel=0;
    int i=0;
    
    for (i=0;i<num_array;i++)
    {
        if ((*(r+i)) >= (r00*gamma_infinity))
        {
            *(gamma+i)=gamma_infinity;
            *(pres+i)=(lumi*pow(r00, 2.0/3.0)*pow(*(r+i), -8.0/3.0) )/(12.0*M_PI*C_LIGHT*pow(gamma_infinity, 4.0/3.0)*pow(C_LIGHT, 2.0)); 
        }
        else
        {
            *(gamma+i)=(*(r+i))/r00;
            *(pres+i)=(lumi*pow(r00, 2.0))/(12.0*M_PI*C_LIGHT*pow(C_LIGHT, 2.0)*pow(*(r+i), 4.0) );  
        }
        
        vel=pow(1-(pow(*(gamma+i), -2.0)) ,0.5);
        *(vx+i)=(vel*(*(x+i)))/pow(pow(*(x+i), 2)+ pow(*(y+i), 2) ,0.5);
        *(vy+i)=(vel*(*(y+i)))/pow(pow(*(x+i), 2)+ pow(*(y+i), 2) ,0.5);
        *(dens+i)=lumi/(4*M_PI*pow(*(r+i), 2.0)*pow(C_LIGHT, 3.0)*gamma_infinity*(*(gamma+i)));
        *(dens_lab+i)=(*(dens+i))*(*(gamma+i));
        *(temp+i)=pow(3*(*(pres+i))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
        //fprintf(fPtr,"Gamma: %lf\nR: %lf\nPres: %e\nvel %lf\nX: %lf\nY %lf\nVx: %lf\nVy: %lf\nDens: %e\nLab_Dens: %e\nTemp: %lf\n", *(gamma+i), *(r+i), *(pres+i), vel, *(x+i), *(y+i), *(vx+i), *(vy+i), *(dens+i), *(dens_lab+i), *(temp+i));

    }
    
}

void structuredFireballPrep(double *r, double *theta,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr)
{
    //This model is provided by Lundman, Peer, Ryde 2014, use this to compare our MCRaT polarization to their polarizations
    double  gamma_0=100, lumi=1e52, r00=1e8, theta_j=1e-2, p=4; //theta_j in paper is 1e-2, 3e-2, 1e-1 and p is 1,2,4
    double T_0=pow(lumi/(4*M_PI*r00*r00*A_RAD*C_LIGHT), 1.0/4.0);
    double eta=0, r_sat=0;
    double vel=0, theta_ratio=0;
    int i=0;
    
    for (i=0;i<num_array;i++)
    {
        
        
        theta_ratio=(*(theta+i))/theta_j;
        eta=gamma_0*pow(1+pow(theta_ratio, 2*p) , -0.5);
        
        if (*(theta+i) >= theta_j*pow(gamma_0/2, 1.0/p))
        {
            //*(gamma+i)=2; //outside with of shear layer have gamma be 2 like in paper
            eta=2.0;
        }
        
        r_sat=eta*r00;
        
        if ((*(r+i)) >= r_sat)
        {
            *(gamma+i)=eta;
            *(temp+i)=T_0*pow(r_sat/(*(r+i)), 2.0/3.0)/eta;
        }
        else
        {
            *(gamma+i)=(*(r+i))/r_sat; //not sure if this is right but it shouldn't matter since we're injecting our photons far from r00
            *(temp+i)=T_0;
        }
        
        vel=pow(1-(pow(*(gamma+i), -2.0)) ,0.5);
        *(vx+i)=(vel*(*(x+i)))/pow(pow(*(x+i), 2)+ pow(*(y+i), 2) ,0.5);
        *(vy+i)=(vel*(*(y+i)))/pow(pow(*(x+i), 2)+ pow(*(y+i), 2) ,0.5);
        *(dens+i)=M_P*lumi/(4*M_PI*M_P*C_LIGHT*C_LIGHT*C_LIGHT*eta*vel*(*(gamma+i))*(*(r+i))*(*(r+i))); //equation paper has extra c, but then units dont work out
        *(dens_lab+i)=(*(dens+i))*(*(gamma+i));
        *(pres+i)=(A_RAD*pow(*(temp+i), 4.0))/(3*pow(C_LIGHT, 2.0));
        //fprintf(fPtr,"eta: %lf\nr_sat: %lf\nGamma: %lf\nR: %lf\nTheta: %lf\nPres: %e\nvel %lf\nX: %lf\nY %lf\nVx: %lf\nVy: %lf\nDens: %e\nLab_Dens: %e\nTemp: %lf\n\n", eta, r_sat, *(gamma+i), *(r+i), (*(theta+i)), *(pres+i), vel, *(x+i), *(y+i), *(vx+i), *(vy+i), *(dens+i), *(dens_lab+i), *(temp+i));
        
    }
    
}


void dirFileMerge(char dir[200], int start_frame, int last_frame, int numprocs, int angle_id, FILE *fPtr )
{
    //function to merge files in mcdir produced by various threads
    double *p0=NULL, *p1=NULL, *p2=NULL, *p3=NULL, *comv_p0=NULL, *comv_p1=NULL, *comv_p2=NULL, *comv_p3=NULL, *r0=NULL, *r1=NULL, *r2=NULL, *s0=NULL, *s1=NULL, *s2=NULL, *s3=NULL, *num_scatt=NULL, *weight=NULL;
    int i=0, j=0, k=0, isNotCorrupted=0, num_types=12;
    int increment=1;
    char filename_k[2000]="", file_no_thread_num[2000]="", cmd[2000]="", mcdata_type[20]="";
    char group[200]="";
    hid_t  file, file_new, group_id, dspace;
    hsize_t dims[1]={0};
    herr_t status, status_group;
    hid_t dset_p0, dset_p1, dset_p2, dset_p3, dset_comv_p0, dset_comv_p1, dset_comv_p2, dset_comv_p3, dset_r0, dset_r1, dset_r2, dset_s0, dset_s1, dset_s2, dset_s3, dset_num_scatt, dset_weight, dset_weight_frame;
   
    //printf("Merging files in %s\n", dir); 
    //#pragma omp parallel for num_threads(num_thread) firstprivate( filename_k, file_no_thread_num, cmd,mcdata_type,num_files, increment ) private(i,j,k)
    // i < last frame because calculation before this function gives last_frame as the first frame of the next process set of frames to merge files for
    
    //if ((COMV_SWITCH!=0) && (STOKES_SWITCH!=0))
    #if COMV_SWITCH == ON && STOKES_SWITCH == ON
    {
        num_types=16;//both switches on, want to save comv and stokes
    }
    //else if ((COMV_SWITCH!=0) || (STOKES_SWITCH!=0))
    #elif COMV_SWITCH == ON || STOKES_SWITCH == ON
    {
        num_types=12;//either switch acivated, just subtract 4 datasets
    }
    #else
    {
        num_types=8;//just save lab 4 momentum, position and num_scatt
    }
    #endif
    
    
    
    for (i=start_frame;i<last_frame;i=i+increment)
    {
        fprintf(fPtr, "Merging files for frame: %d\n", i);
        fflush(fPtr);
        
        #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
        if (i>=3000)
        {
            increment=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
        }
        #endif
        
        j=0;
        for (k=0;k<numprocs;k++)
        {
            //for each process' file, find out how many elements and add up to find total number of elements needed in the data set for the frame number
            snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dir,"mc_proc_", k, ".h5" );
            
            //open the file
            file=H5Fopen(filename_k, H5F_ACC_RDONLY, H5P_DEFAULT);
            
            //see if the frame exists
            snprintf(group,sizeof(group),"%d",i );
            status = H5Eset_auto(NULL, NULL, NULL);
            status_group = H5Gget_objinfo (file, group, 0, NULL);
            status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
            
            //if it does open it and read in the size
            if (status_group == 0)
            {
                //open the datatset
                group_id = H5Gopen2(file, group, H5P_DEFAULT);
                dset_p0 = H5Dopen (group_id, "P0", H5P_DEFAULT); //open dataset
                
                //get the number of points
                dspace = H5Dget_space (dset_p0);
                status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                j+=dims[0];//calculate the total number of photons to save to new hdf5 file
                
                status = H5Sclose (dspace);
                status = H5Dclose (dset_p0);
                status = H5Gclose(group_id);
            }
            status = H5Fclose(file);
        }
        
        //for continuing if the simulation gets stopped, check to see if the new file exists and if the information is correct
        //if the information is incorrect, create file by overwriting it, otherwise dont need to do anything
        snprintf(file_no_thread_num,sizeof(file_no_thread_num),"%s%s%d%s",dir,"mcdata_", i, ".h5" );
        status = H5Eset_auto(NULL, NULL, NULL); //turn off automatic error printing
        file_new=H5Fcreate(file_no_thread_num, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT); //see if the file initially does/doesnt exist
        
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr); //turn on auto error printing

    
        if (file_new<0)
        {
            //fprintf(fPtr, "Checking File %s\n",file_no_thread_num );
            //fflush(fPtr);
            //the file exists, open it with read write 
            file_new=H5Fopen(file_no_thread_num, H5F_ACC_RDWR, H5P_DEFAULT);
            
            for (k=0;k<num_types;k++)
            {
                //if ((COMV_SWITCH!=0) && (STOKES_SWITCH!=0))
                #if COMV_SWITCH == ON && STOKES_SWITCH == ON
                {
                    switch (k)
                    {
                        case 0: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P0"); break;
                        case 1: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P1");break;
                        case 2: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P2"); break;
                        case 3: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P3"); break;
                        case 4: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P0"); break;
                        case 5: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P1");break;
                        case 6: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P2"); break;
                        case 7: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P3"); break;
                        case 8: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R0"); break;
                        case 9: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R1"); break;
                        case 10: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R2"); break;
                        case 11: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S0"); break;
                        case 12: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S1");break;
                        case 13: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S2"); break;
                        case 14: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S3"); break;
                        case 15: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "NS"); break;
                    }
                }
                //else if (STOKES_SWITCH!=0)
                #elif STOKES_SWITCH == ON && COMV_SWITCH == OFF
                {
                    switch (k)
                    {
                        case 0: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P0"); break;
                        case 1: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P1");break;
                        case 2: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P2"); break;
                        case 3: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P3"); break;
                        case 4: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R0"); break;
                        case 5: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R1"); break;
                        case 6: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R2"); break;
                        case 7: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S0"); break;
                        case 8: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S1");break;
                        case 9: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S2"); break;
                        case 10: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "S3"); break;
                        case 11: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "NS"); break;
                    }
                }
                //else if (COMV_SWITCH!=0)
                #elif STOKES_SWITCH == OFF && COMV_SWITCH == ON
                {
                    switch (k)
                    {
                        case 0: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P0"); break;
                        case 1: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P1");break;
                        case 2: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P2"); break;
                        case 3: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P3"); break;
                        case 4: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P0"); break;
                        case 5: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P1");break;
                        case 6: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P2"); break;
                        case 7: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "COMV_P3"); break;
                        case 8: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R0"); break;
                        case 9: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R1"); break;
                        case 10: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R2"); break;
                        case 11: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "NS"); break;
                    }
                }
                #else
                {
                    switch (k)
                    {
                        case 0: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P0"); break;
                        case 1: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P1");break;
                        case 2: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P2"); break;
                        case 3: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "P3"); break;
                        case 4: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R0"); break;
                        case 5: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R1"); break;
                        case 6: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "R2"); break;
                        case 7: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "NS"); break;
                    }
                }
                #endif
            
                //open the datatset
                dset_p0 = H5Dopen (file_new, mcdata_type, H5P_DEFAULT); //open dataset
                
                //get the number of points
                dspace = H5Dget_space (dset_p0);
                status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                
                //fprintf(fPtr, "j:%d, dim: %d\n",j, dims[0] );
                //fflush(fPtr);
                
                isNotCorrupted += fmod(dims[0], j); //if the dimension is the dame then the fmod ==0 (remainder of 0), if all datatsets are ==0 then you get a truth value of 0 meaning that it isnt corrupted
                
                status = H5Sclose (dspace);
                status = H5Dclose (dset_p0);
            }
            
            status = H5Fclose(file_new);
            file_new=-1; //do this so if the file exists it doesnt go into the rewriting portion if the file does exist
        }
        
        //fprintf(fPtr, "file %s has isNotCorrupted=%d\n", file_no_thread_num, isNotCorrupted );
        //fflush(fPtr);
        
        //if the new file doesnt have the dimensions that it should, open it and write over the file, or if the file doesnt exist
        if ((file_new>=0) || (isNotCorrupted != 0 ))
        {
            
            //fprintf(fPtr, "In IF\n" );
            //fflush(fPtr);
            
            if (isNotCorrupted != 0)
            {
                //if the data is corrupted overwrite the file
                file_new = H5Fcreate (file_no_thread_num, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            }
        
            //now allocate enough ememory for j number of points
            p0=malloc(j*sizeof(double));  p1=malloc(j*sizeof(double));  p2=malloc(j*sizeof(double));  p3=malloc(j*sizeof(double));
            comv_p0=malloc(j*sizeof(double));  comv_p1=malloc(j*sizeof(double));  comv_p2=malloc(j*sizeof(double));  comv_p3=malloc(j*sizeof(double));
            r0=malloc(j*sizeof(double));  r1=malloc(j*sizeof(double));  r2=malloc(j*sizeof(double));
            s0=malloc(j*sizeof(double));  s1=malloc(j*sizeof(double));  s2=malloc(j*sizeof(double));  s3=malloc(j*sizeof(double));
            num_scatt=malloc(j*sizeof(double));
            weight=malloc(j*sizeof(double));
        
            j=0;
            for (k=0;k<numprocs;k++)
            {
                //for each process open and read the contents of the dataset
                snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dir,"mc_proc_", k, ".h5" );
                file=H5Fopen(filename_k, H5F_ACC_RDONLY, H5P_DEFAULT);
            
                snprintf(group,sizeof(group),"%d",i );
                status = H5Eset_auto(NULL, NULL, NULL);
                status_group = H5Gget_objinfo (file, group, 0, NULL);
                status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
            
                if (status_group == 0)
                {
                    //open the datatset
                    group_id = H5Gopen2(file, group, H5P_DEFAULT);
                    dset_p0 = H5Dopen (group_id, "P0", H5P_DEFAULT); //open dataset
                    dset_p1 = H5Dopen (group_id, "P1", H5P_DEFAULT);
                    dset_p2 = H5Dopen (group_id, "P2", H5P_DEFAULT);
                    dset_p3 = H5Dopen (group_id, "P3", H5P_DEFAULT);
                    
                    //if (COMV_SWITCH!=0)
                    #if COMV_SWITCH == ON
                    {
                        dset_comv_p0 = H5Dopen (group_id, "COMV_P0", H5P_DEFAULT); //open dataset
                        dset_comv_p1 = H5Dopen (group_id, "COMV_P1", H5P_DEFAULT);
                        dset_comv_p2 = H5Dopen (group_id, "COMV_P2", H5P_DEFAULT);
                        dset_comv_p3 = H5Dopen (group_id, "COMV_P3", H5P_DEFAULT);
                    }
                    #endif
                    
                    dset_r0 = H5Dopen (group_id, "R0", H5P_DEFAULT); 
                    dset_r1 = H5Dopen (group_id, "R1", H5P_DEFAULT);
                    dset_r2 = H5Dopen (group_id, "R2", H5P_DEFAULT);
                    
                    //if (STOKES_SWITCH!=0)
                    #if STOKES_SWITCH == ON
                    {
                        dset_s0 = H5Dopen (group_id, "S0", H5P_DEFAULT);
                        dset_s1 = H5Dopen (group_id, "S1", H5P_DEFAULT);
                        dset_s2 = H5Dopen (group_id, "S2", H5P_DEFAULT);
                        dset_s3 = H5Dopen (group_id, "S3", H5P_DEFAULT);
                    }
                    #endif
                    
                    dset_num_scatt = H5Dopen (group_id, "NS", H5P_DEFAULT);
                    dset_weight = H5Dopen (group_id, "Weight", H5P_DEFAULT); // have to account for this only being used for synchrotron emission switch being on
                
                    //read the data in
                    status = H5Dread(dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p0+j));
                    status = H5Dread(dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p1+j));
                    status = H5Dread(dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p2+j));
                    status = H5Dread(dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p3+j));
                    
                    //if (COMV_SWITCH!=0)
                    #if COMV_SWITCH == ON
                    {
                        status = H5Dread(dset_comv_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (comv_p0+j));
                        status = H5Dread(dset_comv_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (comv_p1+j));
                        status = H5Dread(dset_comv_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (comv_p2+j));
                        status = H5Dread(dset_comv_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (comv_p3+j));
                    }
                    #endif
                    
                    status = H5Dread(dset_r0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (r0+j));
                    status = H5Dread(dset_r1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (r1+j));
                    status = H5Dread(dset_r2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (r2+j));
                    
                    //if (STOKES_SWITCH!=0)
                    #if STOKES_SWITCH == ON
                    {
                        status = H5Dread(dset_s0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (s0+j));
                        status = H5Dread(dset_s1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (s1+j));
                        status = H5Dread(dset_s2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (s2+j));
                        status = H5Dread(dset_s3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (s3+j));
                    }
                    #endif
                    
                    status = H5Dread(dset_num_scatt, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (num_scatt+j));
                    status = H5Dread(dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (weight+j));
                
                    //get the number of points
                    dspace = H5Dget_space (dset_p0);
                    status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                    j+=dims[0];//calculate the total number of photons to save to new hdf5 file
                
                
                    status = H5Sclose (dspace);
                    status = H5Dclose (dset_p0); status = H5Dclose (dset_p1); status = H5Dclose (dset_p2); status = H5Dclose (dset_p3);
                    //if (COMV_SWITCH!=0)
                    #if COMV_SWITCH == ON
                    {
                        status = H5Dclose (dset_comv_p0); status = H5Dclose (dset_comv_p1); status = H5Dclose (dset_comv_p2); status = H5Dclose (dset_comv_p3);
                    }
                    #endif
                    
                    status = H5Dclose (dset_r0); status = H5Dclose (dset_r1); status = H5Dclose (dset_r2);
                    
                    //if (STOKES_SWITCH!=0)
                    #if STOKES_SWITCH == ON
                    {
                        status = H5Dclose (dset_s0); status = H5Dclose (dset_s1); status = H5Dclose (dset_s2); status = H5Dclose (dset_s3);
                    }
                    #endif
                    
                    status = H5Dclose (dset_num_scatt); status = H5Dclose (dset_weight);

                    status = H5Gclose(group_id);
                }
                status = H5Fclose(file);
            }
        
            //create the datatspace and dataset
            dims[0]=j;
            dspace = H5Screate_simple(1, dims, NULL);
            dset_p0=H5Dcreate2(file_new, "P0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_p1=H5Dcreate2(file_new, "P1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_p2=H5Dcreate2(file_new, "P2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_p3=H5Dcreate2(file_new, "P3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            //if (COMV_SWITCH!=0)
            #if COMV_SWITCH == ON
            {
                dset_comv_p0=H5Dcreate2(file_new, "COMV_P0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_comv_p1=H5Dcreate2(file_new, "COMV_P1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_comv_p2=H5Dcreate2(file_new, "COMV_P2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_comv_p3=H5Dcreate2(file_new, "COMV_P3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }
            #endif
            dset_r0=H5Dcreate2(file_new, "R0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_r1=H5Dcreate2(file_new, "R1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_r2=H5Dcreate2(file_new, "R2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            //if (STOKES_SWITCH!=0)
            #if STOKES_SWITCH == ON
            {
                dset_s0=H5Dcreate2(file_new, "S0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_s1=H5Dcreate2(file_new, "S1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_s2=H5Dcreate2(file_new, "S2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_s3=H5Dcreate2(file_new, "S3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }
            #endif
            dset_num_scatt=H5Dcreate2(file_new, "NS", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_weight=H5Dcreate2(file_new, "Weight", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            
            //save the data in the new file
            status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p0);
        
            status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p1);
                        
            status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p2);
                        
            status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p3);
            
            //if (COMV_SWITCH!=0)
            #if COMV_SWITCH == ON
            {
                status = H5Dwrite (dset_comv_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                   H5P_DEFAULT, comv_p0);
                
                status = H5Dwrite (dset_comv_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                   H5P_DEFAULT, comv_p1);
                
                status = H5Dwrite (dset_comv_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                   H5P_DEFAULT, comv_p2);
                
                status = H5Dwrite (dset_comv_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                   H5P_DEFAULT, comv_p3);
            }
            #endif
                        
            status = H5Dwrite (dset_r0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, r0);
        
            status = H5Dwrite (dset_r1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, r1);
                        
            status = H5Dwrite (dset_r2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, r2);
            
            //if (STOKES_SWITCH!=0)
            #if STOKES_SWITCH == ON
            {
                status = H5Dwrite (dset_s0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, s0);
            
                status = H5Dwrite (dset_s1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, s1);
                
                status = H5Dwrite (dset_s2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, s2);
                
                status = H5Dwrite (dset_s3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, s3);
            }
            #endif
                        
            status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, num_scatt);
            
            status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, weight);
        
            status = H5Sclose (dspace);
            status = H5Dclose (dset_p0); status = H5Dclose (dset_p1); status = H5Dclose (dset_p2); status = H5Dclose (dset_p3);
            //if (COMV_SWITCH!=0)
            #if COMV_SWITCH == ON
            {
                status = H5Dclose (dset_comv_p0); status = H5Dclose (dset_comv_p1); status = H5Dclose (dset_comv_p2); status = H5Dclose (dset_comv_p3);
            }
            #endif
            status = H5Dclose (dset_r0); status = H5Dclose (dset_r1); status = H5Dclose (dset_r2);
            //if (STOKES_SWITCH!=0)
            #if STOKES_SWITCH == ON
            {
                status = H5Dclose (dset_s0); status = H5Dclose (dset_s1); status = H5Dclose (dset_s2); status = H5Dclose (dset_s3);
            }
            #endif
            status = H5Dclose (dset_num_scatt); status = H5Dclose (dset_weight);

            status = H5Fclose (file_new);
        
            free(p0);free(p1); free(p2);free(p3);
            free(comv_p0);free(comv_p1); free(comv_p2);free(comv_p3);
            free(r0);free(r1); free(r2);
            free(s0);free(s1); free(s2);free(s3);
            free(num_scatt); free(weight);
        
            isNotCorrupted=0;
        }
        
        
    }
        
    
    
    if (angle_id==0)
    {
        //merge photon weight files
        snprintf(file_no_thread_num,sizeof(file_no_thread_num),"%s%s",dir,"mcdata_PW.h5" );
        file_new = H5Fcreate (file_no_thread_num, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //just recreate the file b/c its easier
        
        //calculate total number of photons 
        j=0;
        for (k=0;k<numprocs;k++)
        {
            //for each process' file, find out how many elements and add up to find total number of elements needed in the data set for the frame number
            snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dir,"mc_proc_", k, ".h5" );
            
            //open the file
            file=H5Fopen(filename_k, H5F_ACC_RDONLY, H5P_DEFAULT);
            
            //open the datatset
            dset_weight = H5Dopen (file, "Weight", H5P_DEFAULT); //open dataset
                
            //get the number of points
            dspace = H5Dget_space (dset_weight);
            status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
            j+=dims[0];//calculate the total number of photons to save to new hdf5 file
                
            status = H5Sclose (dspace);
            status = H5Dclose (dset_weight);
        }
        status = H5Fclose(file);
        
        weight=malloc(j*sizeof(double));
        j=0;
        //open each data set and save it
        for (k=0;k<numprocs;k++)
        {
            snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dir,"mc_proc_", k, ".h5" );
            file=H5Fopen(filename_k, H5F_ACC_RDONLY, H5P_DEFAULT);
            dset_weight = H5Dopen (file, "Weight", H5P_DEFAULT);
            status = H5Dread(dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (weight+j));
            
            dspace = H5Dget_space (dset_weight);
            status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
            j+=dims[0];//calculate the total number of photons to save to new hdf5 file
            
            status = H5Sclose (dspace);
            status = H5Dclose (dset_weight);
            status = H5Fclose(file);
        }
        
        //save the dataset to the new file
        dims[0]=j;
        dspace = H5Screate_simple(1, dims, NULL);
        dset_weight=H5Dcreate2(file_new, "Weight", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, weight);
        
        status = H5Sclose (dspace);
        status = H5Dclose (dset_weight);
        status = H5Fclose(file_new);
        
        free(weight);
    }
    
    //exit(0);
    
}

void modifyFlashName(char flash_file[200], char prefix[200], int frame)
{
    int lim1=0, lim2=0, lim3=0;
    char test[200]="" ;
    //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
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
        lim1=100;
        lim2=1000;
        lim3=10000;
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


void readHydro2D(char hydro_prefix[200], int frame, double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r, double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, FILE *fPtr)
{
    
    FILE *hydroPtr=NULL;
    char hydrofile[200]="", file_num[200]="", full_file[200]="", file_end[200]=""  ;
    char buf[10]="";
    int i=0, j=0, k=0, elem=0, elem_factor=0;
    int all_index_buffer=0, r_min_index=0, r_max_index=0, theta_min_index=0, theta_max_index=0; //all_index_buffer contains phi_min, phi_max, theta_min, theta_max, r_min, r_max indexes to get from grid files
    int r_index=0, theta_index=0;
    float buffer=0;
    float *dens_unprc=NULL,*vel_r_unprc=NULL, *vel_theta_unprc=NULL,*pres_unprc=NULL;
    double ph_rmin=0, ph_rmax=0;
    double r_in=1e10;
    //double *r_edge=malloc(sizeof(double)*(R_DIM_2D+1));
    //double *dr=malloc(sizeof(double)*(R_DIM_2D));
    double *r_unprc=malloc(sizeof(double)*R_DIM_2D);
    double *theta_unprc=malloc(sizeof(double)*THETA_DIM_2D);
    
    if (ph_inj_switch==0)
    {
        ph_rmin=min_r;
        ph_rmax=max_r;
    }
    
    snprintf(file_end,sizeof(file_end),"%s","small.data" );

    //density
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 1,"-" );
    modifyFlashName(file_num, hydrofile, frame);
    
    fprintf(fPtr,">> Opening file %s\n", file_num);
    fflush(fPtr);
    
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&theta_min_index, sizeof(int)*1, 1,hydroPtr);
    fread(&theta_max_index, sizeof(int)*1, 1,hydroPtr);
    fread(&r_min_index, sizeof(int)*1, 1,hydroPtr);
    fread(&r_max_index, sizeof(int)*1, 1,hydroPtr);
    fclose(hydroPtr);
    
    //fortran indexing starts @ 1, but C starts @ 0
    r_min_index--;//=r_min_index-1;
    r_max_index--;//=r_max_index-1;
    theta_min_index--;//=theta_min_index-1;
    theta_max_index--;//=theta_max_index-1;
    
    elem=(r_max_index+1-r_min_index)*(theta_max_index+1-theta_min_index); //max index is max number of elements minus 1, there add one to get total number of elements
    fprintf(fPtr,"Elem %d\n", elem);
    fprintf(fPtr,"Limits %d, %d, %d, %d, %d, %d\n", all_index_buffer, all_index_buffer, theta_min_index, theta_max_index, r_min_index, r_max_index); 
    fflush(fPtr);
    
    //now with number of elements allocate data
    dens_unprc=malloc(elem*sizeof(float));
    vel_r_unprc=malloc(elem*sizeof(float));
    vel_theta_unprc=malloc(elem*sizeof(float));
    pres_unprc=malloc(elem*sizeof(float));
    
    
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    /*
    for (i=0;i<elem;i++)
    {
        fread((dens_unprc+i), sizeof(float),1, hydroPtr); //read  data
    }
    */
    fread(dens_unprc, sizeof(float),elem, hydroPtr);
    fclose(hydroPtr);
    
    //V_r
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 2,"-" );
    modifyFlashName(file_num, hydrofile, frame);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(vel_r_unprc, sizeof(float),elem, hydroPtr); //data
    
    fclose(hydroPtr);
    
    //V_theta
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 3,"-" );
    modifyFlashName(file_num, hydrofile, frame);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    fread(vel_theta_unprc, sizeof(float), elem, hydroPtr); //data
    fclose(hydroPtr);
    
    //u04 is phi component but is all 0
    
    //pres
    snprintf(hydrofile,sizeof(hydrofile),"%s%s%d%s",hydro_prefix,"u0", 8,"-" );
    modifyFlashName(file_num, hydrofile, frame);
    snprintf(full_file, sizeof(full_file), "%s%s", file_num, file_end);
    
    //fprintf(fPtr,">> Opening file %s\n", full_file);
    //fflush(fPtr);
    
    hydroPtr=fopen(full_file, "rb");
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr); //min and max indexes for the grid
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&all_index_buffer, sizeof(int)*1, 1,hydroPtr);
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    fread(&buffer, sizeof(float), 1,hydroPtr); //random stuff about the file from fortran
    
    
    //elem=(r_max_index-r_min_index)*(theta_max_index-theta_min_index);
    //fprintf(fPtr,"Elem %d\n", elem);
    //fprintf(fPtr,"Limits %d, %d, %d, %d, %d, %d\n", all_index_buffer, all_index_buffer, theta_min_index, theta_max_index, r_min_index, r_max_index); 
    //fflush(fPtr);

    fread(pres_unprc, sizeof(float),elem, hydroPtr); //data
   
    fclose(hydroPtr);
    
    
      /*
     for (j=0 ;j<(theta_max_index+1-theta_min_index); j++)
     { 
         for (k=0; k<(r_max_index+1-r_min_index); k++)
         {
         
             fprintf(fPtr,"Pres %d: %e\n", ( j*(r_max_index+1-r_min_index)+k ), *(pres_unprc+( j*(r_max_index+1-r_min_index)+k )));
             //fprintf(fPtr,"Pres %d: %e\n", ( j*(r_max_index)+k ), *(pres_unprc+( j*(r_max_index)+k )));
             fflush(fPtr);
         
         }
     }
         exit(0); 
    */
     
    //R
    snprintf(hydrofile,sizeof(hydrofile),"%s%s",hydro_prefix,"grid-x1.data" );
    hydroPtr=fopen(hydrofile, "r");
    //fprintf(fPtr,">> Opening file %s\n", hydrofile);
    //fflush(fPtr);
    
    i=0;
    while (i<R_DIM_2D)
    {
        fscanf(hydroPtr, "%lf", (r_unprc+i));  //read value
        fgets(buf, 3,hydroPtr); //read comma
        /*
         if (i<5)
         {
             //printf("Here\n");
             fprintf(fPtr,"R %d: %e\n", i, *(r_unprc+i));
             fflush(fPtr);
         }
        */
        i++;
    }
    fclose(hydroPtr);

    
    //theta from y axis
    snprintf(hydrofile,sizeof(hydrofile),"%s%s",hydro_prefix,"grid-x2.data" );
    hydroPtr=fopen(hydrofile, "r");
    //fprintf(fPtr,">> Opening file %s\n", hydrofile);
    //fflush(fPtr);
    
    i=0;
    while (i<THETA_DIM_2D)
    {
        fscanf(hydroPtr, "%lf", (theta_unprc+i));  //read value
        fgets(buf, 3,hydroPtr); //read comma
        /*
        if (i<5)
        {
            fprintf(fPtr,"Theta %d: %e\n", i, *(theta_unprc+i));
            fflush(fPtr);
        }
        */
        i++;
    }
    fclose(hydroPtr);
    
    
    //limit number of array elements
    //fill in radius array and find in how many places r > injection radius
    elem_factor=0;
    elem=0;
    while (elem==0)
    {
        elem_factor++;
        elem=0;
        for (j=0 ;j<(theta_max_index+1-theta_min_index); j++)
        {   
            for (k=0; k<(r_max_index+1-r_min_index); k++)
            {
                i=r_min_index+k; //look at indexes of r that are included in small hydro file
                //if I have photons do selection differently than if injecting photons
                if (ph_inj_switch==0)
                {
                    //if calling this function when propagating photons, choose blocks based on where the photons are
                    if (((ph_rmin - elem_factor*C_LIGHT/fps)<(*(r_unprc+i))) && (*(r_unprc+i)  < (ph_rmax + elem_factor*C_LIGHT/fps) ))
                    {
                        // *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )
                        elem++;
                    }
                }
                else
                {
                    //if calling this function to inject photons choose blocks based on injection parameters, r_inj, which is sufficient
                    if (((r_inj - C_LIGHT/fps)<(*(r_unprc+i))) && (*(r_unprc+i)  < (r_inj + C_LIGHT/fps) ))
                    {
                        // *(pres_unprc+(i*R_DIM*THETA_DIM + j*R_DIM + k  )
                        elem++;
                    }
                
                }
            
            }
        }
    }
    fprintf(fPtr, "Number of post restricted Elems: %d %e\n", elem, r_inj);
    fflush(fPtr);

   
    (*pres)=malloc (elem * sizeof (double ));
    (*velx)=malloc (elem * sizeof (double ));
    (*vely)=malloc (elem * sizeof (double ));
    (*dens)=malloc (elem * sizeof (double ));
    (*x)=malloc (elem * sizeof (double ));
    (*y)=malloc (elem * sizeof (double ));
    (*r)=malloc (elem * sizeof (double ));
    (*theta)=malloc (elem * sizeof (double ));
    (*gamma)=malloc (elem * sizeof (double ));
    (*dens_lab)=malloc (elem * sizeof (double ));
    //szx becomes delta r szy becomes delta theta
    (*szx)=malloc (elem * sizeof (double ));
    (*szy)=malloc (elem * sizeof (double ));
    (*temp)=malloc (elem * sizeof (double ));

    elem=0;
    for (j=0 ;j<(theta_max_index+1-theta_min_index); j++)
     { 
         for (k=0; k<(r_max_index+1-r_min_index); k++)
         {
            r_index=r_min_index+k; //look at indexes of r that are included in small hydro file
            theta_index=theta_min_index+j;
            
            if (ph_inj_switch==0)
            {
                if (((ph_rmin - elem_factor*C_LIGHT/fps)<(*(r_unprc+r_index))) && (*(r_unprc+r_index)  < (ph_rmax + elem_factor*C_LIGHT/fps) ))
                {
                    (*pres)[elem]=*(pres_unprc+( j*(r_max_index+1-r_min_index)+k ));
                    (*velx)[elem]=(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )))*sin(*(theta_unprc+theta_index))+(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )))*cos(*(theta_unprc+theta_index));
                    (*vely)[elem]=(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )))*cos(*(theta_unprc+theta_index))-(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )))*sin(*(theta_unprc+theta_index));
                    (*dens)[elem]=*(dens_unprc+( j*(r_max_index+1-r_min_index)+k ));
                    (*x)[elem]=(*(r_unprc+r_index))*sin(*(theta_unprc+theta_index));
                    (*y)[elem]=(*(r_unprc+r_index))*cos(*(theta_unprc+theta_index));
                    (*r)[elem]=*(r_unprc+r_index);
                    (*szx)[elem]=(*(r_unprc+r_index))*((M_PI/2)/2000);
                    (*szy)[elem]=(M_PI/2)/2000;
                    (*theta)[elem]=*(theta_unprc+theta_index);//theta in radians in relation to jet axis
                    (*gamma)[elem]=pow(pow(1.0-(pow(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )),2)+pow(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )),2)),0.5),-1); //v is in units of c
                    (*dens_lab)[elem]= (*(dens_unprc+( j*(r_max_index+1-r_min_index)+k ))) * pow(pow(1.0-(pow(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )),2)+pow(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )),2)),0.5),-1);
                    (*temp)[elem]=pow(3*(*(pres_unprc+( j*(r_max_index+1-r_min_index)+k )))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                    elem++;
                }
            }
            else
            {
                if (((r_inj - C_LIGHT/fps)<(*(r_unprc+r_index))) && (*(r_unprc+r_index)  < (r_inj + C_LIGHT/fps) ))
                {
                    (*pres)[elem]=*(pres_unprc+( j*(r_max_index+1-r_min_index)+k ));
                    (*velx)[elem]=(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )))*sin(*(theta_unprc+theta_index))+(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )))*cos(*(theta_unprc+theta_index));
                    (*vely)[elem]=(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )))*cos(*(theta_unprc+theta_index))-(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )))*sin(*(theta_unprc+theta_index));
                    (*dens)[elem]=*(dens_unprc+( j*(r_max_index+1-r_min_index)+k ));
                    (*x)[elem]=(*(r_unprc+r_index))*sin(*(theta_unprc+theta_index));
                    (*y)[elem]=(*(r_unprc+r_index))*cos(*(theta_unprc+theta_index));
                    (*r)[elem]=*(r_unprc+r_index);
                    (*szx)[elem]=(*(r_unprc+r_index))*((M_PI/2)/2000);
                    (*szy)[elem]=(M_PI/2)/2000;
                    (*theta)[elem]=*(theta_unprc+theta_index);//theta in radians in relation to jet axis
                    (*gamma)[elem]=pow(pow(1.0-(pow(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )),2)+pow(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )),2)),0.5),-1); //v is in units of c
                    (*dens_lab)[elem]= (*(dens_unprc+( j*(r_max_index+1-r_min_index)+k ))) * pow(pow(1.0-(pow(*(vel_r_unprc+( j*(r_max_index+1-r_min_index)+k )),2)+pow(*(vel_theta_unprc+( j*(r_max_index+1-r_min_index)+k )),2)),0.5),-1);
                    (*temp)[elem]=pow(3*(*(pres_unprc+( j*(r_max_index+1-r_min_index)+k )))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                    elem++;
                    
                }
            }
        }
    }

    (*number)=elem;
    //fprintf(fPtr, "Number of post restricted Elems: %d %e\n", elem, r_inj);
    //fflush(fPtr);
    
    
    free(pres_unprc); //works when not being freed?
    //fprintf(fPtr, "pres Done\n\n");
    //fflush(fPtr);
    
    free(vel_r_unprc);
    //fprintf(fPtr, "vel_r Done\n\n");
    //fflush(fPtr);
    
    free(vel_theta_unprc);
    //fprintf(fPtr, "vel_theta Done\n\n");
    //fflush(fPtr);
    
    free(dens_unprc);
    //fprintf(fPtr, "dens Done\n\n");
    //fflush(fPtr);
    
    free(r_unprc); 
    //fprintf(fPtr, "r Done\n\n");
    //fflush(fPtr);
    
    free(theta_unprc); 
    //fprintf(fPtr, "theta Done\n\n");
    //fflush(fPtr);
    
    pres_unprc=NULL;
    vel_r_unprc=NULL;
    vel_theta_unprc=NULL;
    dens_unprc=NULL;
    r_unprc=NULL;
    theta_unprc=NULL;
    
    //fprintf(fPtr, "ALL Done\n\n");
    //fflush(fPtr);
}

