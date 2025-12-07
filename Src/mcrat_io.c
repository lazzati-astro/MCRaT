//
//  mcrat_io.c
//  
//
//  Created by Tyler Parsotan on 7/23/21. test
//

#include "mcrat.h"

int getOrigNumProcesses(int *counted_cont_procs,  int **proc_array, char dir[STR_BUFFER], int angle_rank,  int angle_procs, int last_frame)
{
    int i=0, j=0, val=0, original_num_procs=-1, rand_num=0;
    int frame2=0, framestart=0, scatt_framestart=0, ph_num=0;
    double time=0;
    char mc_chkpt_files[STR_BUFFER]="", restrt=""; //define new variable that wont write over the restrt variable in the main part of the code, when its put into the readCheckpoint function
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
    //char mc_chkpt_files[STR_BUFFER]="";
    
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

void printPhotons(struct photonList *photon_list, int num_ph_abs, int num_cyclosynch_ph_emit, int num_null_ph, int scatt_cyclosynch_num_ph, int frame,int frame_inj, int frame_last, char dir[STR_BUFFER], int angle_rank, FILE *fPtr )
{
    //function to save the photons' positions and 4 momentum
    
     //now using hdf5 file for each process w/ group structure /(weights or Hydro File #)/(p0,p1,p2,p3, r0, r1, r2, s0, s1, s2, or num_scatt)
    
     
     //open the file if it exists and see if the group exists for the given frame, if frame doesnt exist then write datasets for all photons as extendable
     //if the frame does exist then read information from the prewritten data and then add new data to it as extended chunk
     
     
    int i=0, count=0, rank=1, net_num_ph=photon_list->num_photons; //can have more photons absorbed than emitted, weight_net_num_ph=(frame==frame_inj) ? num_ph-num_ph_abs-num_null_ph : scatt_cyclosynch_num_ph
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    char mc_file[STR_BUFFER]="", group[200]="", group_weight[200]="", *ph_type=NULL;
    double p0[net_num_ph], p1[net_num_ph], p2[net_num_ph], p3[net_num_ph] , r0[net_num_ph], r1[net_num_ph], r2[net_num_ph], num_scatt[net_num_ph], weight[net_num_ph], global_weight[net_num_ph];
    double s0[net_num_ph], s1[net_num_ph], s2[net_num_ph], s3[net_num_ph], comv_p0[net_num_ph], comv_p1[net_num_ph], comv_p2[net_num_ph], comv_p3[net_num_ph];
    hid_t  file, file_init, dspace, dspace_weight, dspace_global_weight, fspace, mspace, prop, prop_weight, prop_global_weight, group_id;
    hid_t dset_p0, dset_p1, dset_p2, dset_p3, dset_r0, dset_r1, dset_r2, dset_s0, dset_s1, dset_s2, dset_s3, dset_num_scatt, dset_weight, dset_weight_2, dset_comv_p0, dset_comv_p1, dset_comv_p2, dset_comv_p3, dset_ph_type;
    herr_t status, status_group, status_weight, status_weight_2;
    hsize_t dims[1]={net_num_ph}, dims_weight[1]={net_num_ph}, dims_old[1]={0}; //1 is the number of dimansions for the dataset, called rank
    struct photon *ph=NULL; //pointer to a photon struct


    
    hsize_t maxdims[1]={H5S_UNLIMITED};
    hsize_t      size[1];
    hsize_t      offset[1];
    
    fprintf(fPtr, "num_ph %d num_ph_abs %d num_null_ph %d num_cyclosynch_ph_emit %d\nAllocated weight to be %d values large and other arrays to be %d\n",photon_list->num_photons,num_ph_abs,num_null_ph,num_cyclosynch_ph_emit, net_num_ph, net_num_ph);
    
    ph_type=malloc((net_num_ph)*sizeof(char));
    
    //save photon data into large arrays, NEED TO KNOW HOW MANY NULL PHOTONS WE HAVE AKA SAVED SPACE THAT AREN'T ACTUALLY PHOTONS TO PROPERLY SAVE SPACE FOR ARRAYS ABOVE
    count=0;//used to keep track of weight values since it may not be the same as num_ph
    //#pragma omp parallel for num_threads(num_thread) reduction(+:weight_net_num_ph)
    for (i=0;i<photon_list->list_capacity;i++)
    {
        ph=getPhoton(photon_list, i);
        
        if (ph->weight != 0)
        {
            p0[count]= (ph->p0);
            p1[count]= (ph->p1);
            p2[count]= (ph->p2);
            p3[count]= (ph->p3);
            r0[count]= (ph->r0);
            r1[count]= (ph->r1);
            r2[count]= (ph->r2);
            #if COMV_SWITCH == ON
            {
                comv_p0[count]= (ph->comv_p0);
                comv_p1[count]= (ph->comv_p1);
                comv_p2[count]= (ph->comv_p2);
                comv_p3[count]= (ph->comv_p3);
            }
            #endif
            #if STOKES_SWITCH == ON
            {
                s0[count]= (ph->s0);
                s1[count]= (ph->s1);
                s2[count]= (ph->s2);
                s3[count]= (ph->s3);
            }
            #endif
            num_scatt[count]= (ph->num_scatt);
            weight[count]= (ph->weight);
             //fprintf(fPtr, "%d %c %e %e %e %e %e %e %e %e\n", i, ph->type, ph->r0, ph->r1, ph->r2, ph->num_scatt, ph->weight, ph->p0, ph->comv_p0, ph->p0*C_LIGHT/1.6e-9);
            
            if ((frame==frame_last))
            {
                global_weight[count]=(ph->weight);
            }
            
            *(ph_type+count)=ph->type;
            //printf("%d %c %e %e %e %e %e %e %e %e %c\n", i, ph->type, ph->r0, ph->r1, ph->r2, ph->num_scatt, ph->weight, ph->p0, ph->comv_p0, ph->p0*C_LIGHT/1.6e-9, *(ph_type+count));
            
            count++;
        }
        
    }
    
    
    //make strings for file name and group
    snprintf(mc_file,sizeof(mc_file),"%s%s%d%s",dir,"mc_proc_", angle_rank, ".h5" );
    snprintf(group,sizeof(group),"%d",frame );
    
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
        //snprintf(group_weight,sizeof(group_weight),"/PW",i );
        status = H5Eset_auto(NULL, NULL, NULL);
        status_weight = H5Gget_objinfo (file, "/PW", 0, NULL);
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
        
        fprintf(fPtr,"Status of /PW %d\n", status_weight);
        
        //the file has been newly created or if the group does not exist then  create the group for the frame
        group_id = H5Gcreate2(file, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        
        /* Modify dataset creation properties, i.e. enable chunking  */
        prop = H5Pcreate (H5P_DATASET_CREATE);
        status = H5Pset_chunk (prop, rank, dims);

        prop_weight= H5Pcreate (H5P_DATASET_CREATE);
        status = H5Pset_chunk (prop_weight, rank, dims_weight);

    
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
        
        #if SAVE_TYPE == ON
        {
            dset_ph_type = H5Dcreate2 (group_id, "PT", H5T_NATIVE_CHAR, dspace,
                                       H5P_DEFAULT, prop, H5P_DEFAULT);
        }
        #endif
        
        dset_num_scatt = H5Dcreate2 (group_id, "NS", H5T_NATIVE_DOUBLE, dspace,
                            H5P_DEFAULT, prop, H5P_DEFAULT);
                            
        dset_weight_2 = H5Dcreate2 (group_id, "PW", H5T_NATIVE_DOUBLE, dspace_weight,
                            H5P_DEFAULT, prop_weight, H5P_DEFAULT); //save the new injected photons' weights

                         
        /* Write data to dataset */
        status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p0);
        
        status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p1);
                        
        status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p2);
                        
        status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, p3);
        
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
        
        #if SAVE_TYPE == ON
        {
            status = H5Dwrite (dset_ph_type, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, ph_type);
        }
        #endif

        
        status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, num_scatt);
        
        status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, weight);
        
        status = H5Pclose (prop_weight);
        status = H5Dclose (dset_weight_2);

        
        status = H5Pclose (prop);
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
        
        #if SAVE_TYPE == ON
        {
            dset_ph_type = H5Dopen (group_id, "PT", H5P_DEFAULT); //open dataset
            dspace = H5Dget_space (dset_ph_type);
            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims
            size[0] = dims[0]+ dims_old[0];
            status = H5Dset_extent (dset_ph_type, size);
            fspace = H5Dget_space (dset_ph_type);
            offset[0] = dims_old[0];
            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL,
                                      dims, NULL);
            mspace = H5Screate_simple (rank, dims, NULL);
            status = H5Dwrite (dset_ph_type, H5T_NATIVE_CHAR, mspace, fspace,
                                H5P_DEFAULT, ph_type);
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
        snprintf(group_weight,sizeof(group_weight),"PW",i );
        status = H5Eset_auto(NULL, NULL, NULL);
        status_weight = H5Gget_objinfo (group_id, "PW", 0, NULL);
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
        
        fprintf(fPtr,"Status of /frame/PW %d\n", status_weight);
        
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
        if (status_weight >= 0)
        {
            //will have to create the weight dataset for the new set of phtons that have been injected, although it may already be created since emitting photons now
            //see if the group exists
            status = H5Eset_auto(NULL, NULL, NULL);
            status_weight_2 = H5Gget_objinfo (group_id, "PW", 0, NULL);
            status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
            
            if (status_weight_2 < 0)
            {
                //the dataset doesnt exist
                 /* Modify dataset creation properties, i.e. enable chunking  */
                prop = H5Pcreate (H5P_DATASET_CREATE);
                status = H5Pset_chunk (prop, rank, dims);
        
                /* Create the data space with unlimited dimensions. */
                dspace = H5Screate_simple (rank, dims, maxdims);
                
                dset_weight_2 = H5Dcreate2 (group_id, "PW", H5T_NATIVE_DOUBLE, dspace,
                                H5P_DEFAULT, prop, H5P_DEFAULT); //save the new injected photons' weights
                
                status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, weight);
                
                status = H5Pclose (prop);
            }
            else
            {
                //it exists and need to modify it
                dset_weight_2 = H5Dopen (group_id, "PW", H5P_DEFAULT); //open dataset
                
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
            
            dset_weight_2 = H5Dcreate2 (group_id, "PW", H5T_NATIVE_DOUBLE, dspace_weight,
                                        H5P_DEFAULT, prop_weight, H5P_DEFAULT);
            
            status = H5Dwrite (dset_weight_2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, weight);
            
            status = H5Pclose (prop_weight);
        }
        
        status = H5Dclose (dset_weight_2);
        
        status = H5Sclose (dspace);
        status = H5Sclose (mspace);
        status = H5Sclose (fspace);
        
    }
    
    
    /* Close resources */
    
    free(ph_type);
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
    
    #if SAVE_TYPE == ON
    {
        status = H5Dclose (dset_ph_type);
    }
    #endif
    
    status = H5Dclose (dset_num_scatt);
    
    /* Close the group. */
    status = H5Gclose(group_id);
    
    /* Terminate access to the file. */
    status = H5Fclose(file);
    

}

int saveCheckpoint(char dir[STR_BUFFER], int frame, int frame2, int scatt_frame, double time_now, struct photonList *photon_list, int last_frame, int angle_rank,int angle_size )
{
    //function to save data necessary to restart simulation if it ends
    //need to save all photon data
    FILE *fPtr=NULL;
    char checkptfile[2000]="";
    char command[2000]="";
    char restart;
    int i=0, success=0, ph_num=0;
    struct photon *ph=NULL;
    
    
    snprintf(checkptfile,sizeof(checkptfile),"%s%s%d%s",dir,"mc_chkpt_", angle_rank,".dat" );
    //snprintf(checkptfile,sizeof(checkptfile),"%s%s%d%s%d%s",dir,"mc_chkpt_", angle_rank, "_frame_", scatt_frame, ".dat" ); //look at frame 1341?

    
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
            restart=CONTINUE;
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
            ph_num=photon_list->list_capacity;
            fwrite(&ph_num, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote ph_num\n",  angle_rank);
            fflush(stdout);
            for(i=0;i<photon_list->list_capacity;i++)
            {
                ph=getPhoton(photon_list, i);

                #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((ph->type == COMPTONIZED_PHOTON) && (ph->weight != 0))
                {
                    ph->type = UNABSORBED_CS_PHOTON; //set this to be an old synchrotron scattered photon
                }
                #endif
                fwrite(ph, sizeof(struct photon ), 1, fPtr);
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
            restart=CONTINUE;
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
            ph_num=photon_list->list_capacity;
            fwrite(&ph_num, sizeof(int), 1, fPtr);
            //printf("Rank: %d wrote ph_num\n",  angle_rank);
            fflush(stdout);
            for(i=0;i<photon_list->list_capacity;i++)
            {
                ph=getPhoton(photon_list, i);
                
                #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((ph->type == COMPTONIZED_PHOTON) && (ph->weight != 0))
                {
                    ph->type = UNABSORBED_CS_PHOTON; //set this to be an old synchrotron scattered photon
                }
                #endif
                //fwrite((ph), sizeof(struct photon )*ph_num, ph_num, fPtr);
                fwrite(ph, sizeof(struct photon ), 1, fPtr);
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
            restart=INITALIZE;
            fwrite(&restart, sizeof(char), 1, fPtr);
            fwrite(&frame, sizeof(int), 1, fPtr);
            fwrite(&frame2, sizeof(int), 1, fPtr);
            for(i=0;i<photon_list->list_capacity;i++)
            {
                ph=getPhoton(photon_list, i);
                #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((ph->type == COMPTONIZED_PHOTON) && (ph->weight != 0))
                {
                    ph->type = UNABSORBED_CS_PHOTON; //set this to be an old synchrotron scattered photon
                }
                #endif
                fwrite(ph, sizeof(struct photon ), 1, fPtr);
            }

            success=0;
        }
    }
    if (success==0)
    {
        fclose(fPtr);
    }
    
    return success;
}

int readCheckpoint(char dir[STR_BUFFER], struct photonList *photon_list, int *frame2, int *framestart, int *scatt_framestart, char *restart, double *time, int angle_rank, int *angle_size )
{
    //function to read in data from checkpoint file
    FILE *fPtr=NULL;
    char checkptfile[STR_BUFFER]="";
    int i=0, ph_num=0;
    int scatt_cyclosynch_num_ph=0;//count the number of scattered synchrotron photons from the previosu frame that were saved
    //int frame, scatt_frame, ph_num, i=0;
    struct photon *phHolder=NULL; //pointer to struct to hold data read in from checkpoint file
    struct photon *ph=NULL; //array of photons that are read in
    
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
        
        if((*restart)==CONTINUE)
        {
            fread(scatt_framestart, sizeof(int), 1, fPtr);
            
            //if ((riken_switch==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && ((*scatt_framestart)>=3000))
            #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
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
            fread(&ph_num, sizeof(int), 1, fPtr);
            //printf("%d\n", *ph_num);
            
            phHolder=malloc(sizeof(struct photon));
            (ph)=malloc(sizeof(struct photon)*(ph_num)); //allocate memory to hold photon data
            
            
            for (i=0;i<(*ph_num);i++)
            {
                fread(phHolder, sizeof(struct photon), 1, fPtr);
                //printf("%e,%e,%e, %e,%e,%e, %e, %e\n",(ph)->p0, (ph)->p1, (ph)->p2, ph->p3, (ph)->r0, (ph)->r1, (ph)->r2, ph->num_scatt );
                
                ph[i].p0=phHolder->p0;
                ph[i].p1=phHolder->p1;
                ph[i].p2=phHolder->p2;
                ph[i].p3=phHolder->p3;
                ph[i].comv_p0=phHolder->comv_p0;
                ph[i].comv_p1=phHolder->comv_p1;
                ph[i].comv_p2=phHolder->comv_p2;
                ph[i].comv_p3=phHolder->comv_p3;
                ph[i].r0= phHolder->r0;
                ph[i].r1=phHolder->r1 ;
                ph[i].r2=phHolder->r2;
                ph[i].s0=phHolder->s0;
                ph[i].s1=phHolder->s1;
                ph[i].s2=phHolder->s2;
                ph[i].s3=phHolder->s3;
                ph[i].num_scatt=phHolder->num_scatt;
                ph[i].weight=phHolder->weight;
                ph[i].nearest_block_index= phHolder->nearest_block_index;
                ph[i].type= phHolder->type;
                
                #if CYCLOSYNCHROTRON_SWITCH == ON
                    if ((ph[i].weight != 0) && ((ph[i].type == COMPTONIZED_PHOTON) || (ph[i].type == UNABSORBED_CS_PHOTON)) && (ph[i].p0 > 0))
                    {
                        scatt_cyclosynch_num_ph++;
                    }
                //printf("%d %c %e %e %e %e %e %e %e\n", i, ph[i].type, ph[i].r0, ph[i].r1, ph[i].r2, ph[i].num_scatt, ph[i].weight, ph[i].p0*C_LIGHT/1.6e-9, ph[i].comv_p0);
                #endif
            }
            
            //szve the whole array to our photon list struct
            setPhotonList(photon_list, ph, ph_num);
            
            free(phHolder);
            free(ph);
            //printf("In readcheckpoint count=%d\n", count);
        }
        else
        {
            //if ((riken_switch==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && ((*framestart)>=3000))
            #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
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
        *restart=INITALIZE;
        
    }
    
    return scatt_cyclosynch_num_ph;
}

void readMcPar(struct hydro_dataframe *hydro_data, double *theta_jmin, double *theta_j, double *n_theta_j, double **inj_radius, int **frm0, int **frm2, int *min_photons, int *max_photons, char *spect, char *restart)
{
    //function to read mc.par file
    char mc_file[STR_BUFFER]="" ;
    FILE *fptr=NULL;
    char buf[STR_BUFFER]="", buf2[STR_BUFFER]="", *value, *context = NULL, copied_str[STR_BUFFER]="";
    double theta_deg;
    int i, val;
    
    //open file
    snprintf(mc_file,sizeof(mc_file),"%s%s%s",FILEPATH, MC_PATH,MCPAR);
    printf(">> MCRaT:  Reading parameter file %s\n", mc_file);
    fptr=fopen(mc_file,"r");
    
    //read first block about hydro simulation frame
    fgets(buf, sizeof(buf), fptr); //reads block info
    fgets(buf, sizeof(buf),fptr); //reads /n
    fscanf(fptr, "%lf", &(hydro_data->fps));
    fgets(buf, sizeof(buf),fptr); //reads until end of line
    fscanf(fptr, "%d",&(hydro_data->last_frame));
    fgets(buf, sizeof(buf),fptr); //reads until end of line
    fscanf(fptr, "%lf", &((hydro_data->r0_domain)[0]) );
    fscanf(fptr, "%lf", &((hydro_data->r0_domain)[1]) );
    fgets(buf, sizeof(buf),fptr); //reads until end of line
    fscanf(fptr, "%lf", &((hydro_data->r1_domain)[0]));
    fscanf(fptr, "%lf", &((hydro_data->r1_domain)[1]));
    fgets(buf, sizeof(buf),fptr); //reads until end of line
    fscanf(fptr, "%lf", &((hydro_data->r2_domain)[0]));
    fscanf(fptr, "%lf", &((hydro_data->r2_domain)[1]));
    fgets(buf, sizeof(buf),fptr); //reads until end of line

    //read second block about MCRaT injection angles
    fgets(buf, sizeof(buf),fptr); //reads block info
    fgets(buf, sizeof(buf),fptr); //reads /n
    fscanf(fptr, "%lf",&theta_deg);
    *theta_jmin=theta_deg;// leave as degrees to manipulate processes
    fgets(buf, sizeof(buf),fptr);
    fscanf(fptr, "%lf",&theta_deg);
    *theta_j=theta_deg;//leave as degrees to manipulate processes
    fgets(buf, sizeof(buf),fptr);
    fscanf(fptr, "%lf",&theta_deg);
    *n_theta_j=theta_deg;
    fgets(buf, sizeof(buf),fptr); //reads the rest of the line
    
    //need to read in next line with n_theta_j values for the injection frame start
    (*inj_radius)=malloc( ((int) *n_theta_j)*sizeof(double) );
    (*frm0)=malloc(((int) *n_theta_j)*sizeof(int));
    (*frm2)=malloc(((int) *n_theta_j)*sizeof(int));
    
    fgets(buf, sizeof(buf),fptr); //reads the whole line for injection frame start
    value = strtok_r(buf, " ", &context);
    for (i=0;i< (int) *n_theta_j;i++)
    {
        strcpy(copied_str, value);
        //printf("i %d Read token: %s\n", i, value);
        (*frm0)[i]=strtol(copied_str, buf2, 10);
        value = strtok_r(NULL, " ", &context);
    }
    
    fgets(buf, sizeof(buf),fptr); //reads the whole line for injection frame end
    value = strtok_r(buf, " ", &context);
    for (i=0;i< (int) *n_theta_j;i++)
    {
        strcpy(copied_str, value);
        //printf("i %d Read token: %s\n", i, value);
        (*frm2)[i]=strtol(copied_str, buf2, 10)+(*frm0)[i];
        value = strtok_r(NULL, " ", &context);
    }
    
    fgets(buf, sizeof(buf),fptr); //reads the whole line for injection frame start
    value = strtok_r(buf, " ", &context);
    for (i=0;i< (int) *n_theta_j;i++)
    {
        strcpy(copied_str, value);
        //printf("i %d Read token: %s\n", i, value);
        (*inj_radius)[i]=strtof(copied_str, NULL);
        value = strtok_r(NULL, " ", &context);
    }
    fgets(buf, sizeof(buf),fptr); //reads the new line
    
    //look at the photon block
    fgets(buf, sizeof(buf),fptr); //reads the block header
    fgets(buf, sizeof(buf),fptr); //reads the \n
    *spect=getc(fptr);
    fgets(buf, sizeof(buf),fptr); //reads the remainder of the line
    
    fscanf(fptr, "%d",min_photons);
    fgets(buf, 100,fptr);
    
    fscanf(fptr, "%d",max_photons);
    fgets(buf, 100,fptr);
    fgets(buf, 100,fptr);
    
    //read the initialization or continuation char
    fgets(buf, sizeof(buf),fptr); //reads the block header
    fgets(buf, sizeof(buf),fptr); //reads the \n
    *restart=getc(fptr);
    fgets(buf, 100,fptr);
    
    //close file
    fclose(fptr);
}

void dirFileMerge(char dir[STR_BUFFER], int start_frame, int last_frame, int numprocs, int angle_id, FILE *fPtr )
{
    //function to merge files in mcdir produced by various threads
    double *p0=NULL, *p1=NULL, *p2=NULL, *p3=NULL, *comv_p0=NULL, *comv_p1=NULL, *comv_p2=NULL, *comv_p3=NULL, *r0=NULL, *r1=NULL, *r2=NULL, *s0=NULL, *s1=NULL, *s2=NULL, *s3=NULL, *num_scatt=NULL, *weight=NULL;
    int i=0, j=0, k=0, isNotCorrupted=0, num_types=9; //just save lab 4 momentum, position and num_scatt by default
    int increment=1;
    char filename_k[STR_BUFFER]="", file_no_thread_num[STR_BUFFER]="", cmd[STR_BUFFER]="", mcdata_type[20]="";
    char group[200]="", *ph_type=NULL;
    hid_t  file, file_new, group_id, dspace;
    hsize_t dims[1]={0};
    herr_t status, status_group;
    hid_t dset_p0, dset_p1, dset_p2, dset_p3, dset_comv_p0, dset_comv_p1, dset_comv_p2, dset_comv_p3, dset_r0, dset_r1, dset_r2, dset_s0, dset_s1, dset_s2, dset_s3, dset_num_scatt, dset_weight, dset_weight_frame, dset_ph_type;
   
    //printf("Merging files in %s\n", dir);
    //#pragma omp parallel for num_threads(num_thread) firstprivate( filename_k, file_no_thread_num, cmd,mcdata_type,num_files, increment ) private(i,j,k)
    // i < last frame because calculation before this function gives last_frame as the first frame of the next process set of frames to merge files for
    
    #if COMV_SWITCH == ON && STOKES_SWITCH == ON
    {
        num_types=17;//both switches on, want to save comv and stokes
    }
    #elif COMV_SWITCH == ON || STOKES_SWITCH == ON
    {
        num_types=13;//either switch acivated, just subtract 4 datasets
    }
    #else
    {
        num_types=9;//just save lab 4 momentum, position and num_scatt
    }
    #endif
    
    #if SAVE_TYPE == ON
    {
        num_types+=1;
    }
    #endif
    
    
    for (i=start_frame;i<last_frame;i=i+increment)
    {
        fprintf(fPtr, "Merging files for frame: %d\n", i);
        fflush(fPtr);
        
        #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
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
                        case 16: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PW"); break;
                        #if SAVE_TYPES == ON
                        {
                            case 17: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PT"); break;
                        }
                        #endif
                    }
                }
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
                        case 12: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PW"); break;
                        #if SAVE_TYPES == ON
                        {
                            case 13: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PT"); break;
                        }
                        #endif
                    }
                }
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
                        case 12: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PW"); break;
                        #if SAVE_TYPES == ON
                        {
                            case 13: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PT"); break;
                        }
                        #endif
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
                        case 8: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PW"); break;
                        #if SAVE_TYPES == ON
                        {
                            case 9: snprintf(mcdata_type,sizeof(mcdata_type), "%s", "PT"); break;
                        }
                        #endif
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
            num_scatt=malloc(j*sizeof(double)); weight=malloc(j*sizeof(double));
            ph_type=malloc((j)*sizeof(char));
        
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
                    
                    #if STOKES_SWITCH == ON
                    {
                        dset_s0 = H5Dopen (group_id, "S0", H5P_DEFAULT);
                        dset_s1 = H5Dopen (group_id, "S1", H5P_DEFAULT);
                        dset_s2 = H5Dopen (group_id, "S2", H5P_DEFAULT);
                        dset_s3 = H5Dopen (group_id, "S3", H5P_DEFAULT);
                    }
                    #endif
                    
                    dset_num_scatt = H5Dopen (group_id, "NS", H5P_DEFAULT);
                    
                    dset_weight = H5Dopen (group_id, "PW", H5P_DEFAULT); // have to account for this only being used for synchrotron emission switch being on
                    
                    #if SAVE_TYPE == ON
                    {
                        dset_ph_type = H5Dopen (group_id, "PT", H5P_DEFAULT);
                    }
                    #endif

                    
                    //read the data in
                    status = H5Dread(dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p0+j));
                    status = H5Dread(dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p1+j));
                    status = H5Dread(dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p2+j));
                    status = H5Dread(dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (p3+j));
                    
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
                    
                    #if SAVE_TYPE == ON
                    {
                        status = H5Dread(dset_ph_type, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, (ph_type+j));
                    }
                    #endif

                    
                
                    //get the number of points
                    dspace = H5Dget_space (dset_p0);
                    status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                    j+=dims[0];//calculate the total number of photons to save to new hdf5 file
                
                
                    status = H5Sclose (dspace);
                    status = H5Dclose (dset_p0); status = H5Dclose (dset_p1); status = H5Dclose (dset_p2); status = H5Dclose (dset_p3);

                    #if COMV_SWITCH == ON
                    {
                        status = H5Dclose (dset_comv_p0); status = H5Dclose (dset_comv_p1); status = H5Dclose (dset_comv_p2); status = H5Dclose (dset_comv_p3);
                    }
                    #endif
                    
                    status = H5Dclose (dset_r0); status = H5Dclose (dset_r1); status = H5Dclose (dset_r2);
                    
                    #if STOKES_SWITCH == ON
                    {
                        status = H5Dclose (dset_s0); status = H5Dclose (dset_s1); status = H5Dclose (dset_s2); status = H5Dclose (dset_s3);
                    }
                    #endif
                    
                    #if SAVE_TYPE == ON
                    {
                        status = H5Dclose (dset_ph_type);
                    }
                    #endif

                    
                    status = H5Dclose (dset_num_scatt);
                    
                    status = H5Dclose (dset_weight);

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
            
            #if STOKES_SWITCH == ON
            {
                dset_s0=H5Dcreate2(file_new, "S0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_s1=H5Dcreate2(file_new, "S1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_s2=H5Dcreate2(file_new, "S2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_s3=H5Dcreate2(file_new, "S3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }
            #endif
            
            dset_num_scatt=H5Dcreate2(file_new, "NS", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
            dset_weight=H5Dcreate2(file_new, "PW", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
            #if SAVE_TYPE == ON
            {
                dset_ph_type=H5Dcreate2(file_new, "PT", H5T_NATIVE_CHAR, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }
            #endif

            
            //save the data in the new file
            status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p0);
        
            status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p1);
                        
            status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p2);
                        
            status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, p3);
            
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
            
            #if SAVE_TYPE == ON
            {
                status = H5Dwrite (dset_ph_type, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, ph_type);
            }
            #endif

                        
            status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, num_scatt);
            
                status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, weight);
            
            status = H5Sclose (dspace);
            status = H5Dclose (dset_p0); status = H5Dclose (dset_p1); status = H5Dclose (dset_p2); status = H5Dclose (dset_p3);
            
            #if COMV_SWITCH == ON
            {
                status = H5Dclose (dset_comv_p0); status = H5Dclose (dset_comv_p1); status = H5Dclose (dset_comv_p2); status = H5Dclose (dset_comv_p3);
            }
            #endif
            
            status = H5Dclose (dset_r0); status = H5Dclose (dset_r1); status = H5Dclose (dset_r2);
            
            #if STOKES_SWITCH == ON
            {
                status = H5Dclose (dset_s0); status = H5Dclose (dset_s1); status = H5Dclose (dset_s2); status = H5Dclose (dset_s3);
            }
            #endif
            
            #if SAVE_TYPE == ON
            {
                status = H5Dclose (dset_ph_type);
            }
            #endif

            status = H5Dclose (dset_num_scatt);
            
            status = H5Dclose (dset_weight);
            

            status = H5Fclose (file_new);
        
            free(p0);free(p1); free(p2);free(p3);
            free(comv_p0);free(comv_p1); free(comv_p2);free(comv_p3);
            free(r0);free(r1); free(r2);
            free(s0);free(s1); free(s2);free(s3);
            free(num_scatt); free(weight);
            free(ph_type);
        
            isNotCorrupted=0;
        }
        
        
    }
        
    
    //exit(0);
    
}

void hydroDataFrameInitialize(struct hydro_dataframe *hydro_data)
{
    //initialize pointers in hydro dataframe to NULL for debugging
    hydro_data->r0=NULL;
    hydro_data->r1=NULL;
    hydro_data->r2=NULL;
    hydro_data->r0_size=NULL;
    hydro_data->r1_size=NULL;
    hydro_data->r2_size=NULL;
    hydro_data->r=NULL;
    hydro_data->theta=NULL;
    hydro_data->v0=NULL;
    hydro_data->v1=NULL;
    hydro_data->v2=NULL;
    hydro_data->dens=NULL;
    hydro_data->dens_lab=NULL;
    hydro_data->pres=NULL;
    hydro_data->temp=NULL;
    hydro_data->gamma=NULL;
    hydro_data->B0=NULL;
    hydro_data->B1=NULL;
    hydro_data->B2=NULL;
    #if NONTHERMAL_E_DIST != OFF
        hydro_data->nonthermal_dens=NULL;
    #endif

}

void freeHydroDataFrame(struct hydro_dataframe *hydro_data)
{
    //free pointers in hydro dataframe to NULL for debugging
    free(hydro_data->r0);
    free(hydro_data->r1);
    free(hydro_data->r2);
    free(hydro_data->r0_size);
    free(hydro_data->r1_size);
    free(hydro_data->r2_size);
    free(hydro_data->r);
    free(hydro_data->theta);
    free(hydro_data->v0);
    free(hydro_data->v1);
    free(hydro_data->v2);
    free(hydro_data->dens);
    free(hydro_data->dens_lab);
    free(hydro_data->pres);
    free(hydro_data->temp);
    free(hydro_data->gamma);
    free(hydro_data->B0);
    free(hydro_data->B1);
    free(hydro_data->B2);
    
    hydro_data->r0=NULL;
    hydro_data->r1=NULL;
    hydro_data->r2=NULL;
    hydro_data->r0_size=NULL;
    hydro_data->r1_size=NULL;
    hydro_data->r2_size=NULL;
    hydro_data->r=NULL;
    hydro_data->theta=NULL;
    hydro_data->v0=NULL;
    hydro_data->v1=NULL;
    hydro_data->v2=NULL;
    hydro_data->dens=NULL;
    hydro_data->dens_lab=NULL;
    hydro_data->pres=NULL;
    hydro_data->temp=NULL;
    hydro_data->gamma=NULL;
    hydro_data->B0=NULL;
    hydro_data->B1=NULL;
    hydro_data->B2=NULL;

    #if NONTHERMAL_E_DIST != OFF
        free(hydro_data->nonthermal_dens);
        hydro_data->nonthermal_dens=NULL;
    #endif

}

void allocateHydroDataFrameMemory(struct hydro_dataframe *hydro_data, int n_elements)
{

    hydro_data->r0=malloc (n_elements * sizeof (double ));
    hydro_data->r1=malloc (n_elements * sizeof (double ));
    hydro_data->r0_size=malloc (n_elements * sizeof (double ));
    hydro_data->r1_size=malloc (n_elements * sizeof (double ));
    hydro_data->r=malloc (n_elements * sizeof (double ));
    hydro_data->theta=malloc (n_elements * sizeof (double ));
    hydro_data->v0=malloc (n_elements * sizeof (double ));
    hydro_data->v1=malloc (n_elements * sizeof (double ));
    hydro_data->dens=malloc (n_elements * sizeof (double ));
    hydro_data->dens_lab=malloc (n_elements * sizeof (double ));
    hydro_data->pres=malloc (n_elements * sizeof (double ));
    hydro_data->temp=malloc (n_elements * sizeof (double ));
    hydro_data->gamma=malloc (n_elements * sizeof (double ));

    #if NONTHERMAL_E_DIST != OFF
        hydro_data->nonthermal_dens=malloc (n_elements * sizeof (double ));
    #endif

    #if B_FIELD_CALC == SIMULATION
        (hydro_data->B0)= malloc (n_elements * sizeof (double));
        (hydro_data->B1)= malloc (n_elements * sizeof (double));
    #endif


    #if DIMENSIONS == THREE
        (hydro_data->r2)=malloc(n_elements*sizeof (double));
        (hydro_data->r2_size)=malloc(n_elements*sizeof (double));
    #endif

    #if DIMENSIONS == THREE || DIMENSIONS == TWO_POINT_FIVE
        (hydro_data->v2)=malloc (n_elements * sizeof (double));
        #if B_FIELD_CALC==SIMULATION
            (hydro_data->B2)= malloc (n_elements * sizeof (double));
        #endif
    #endif


}

int getHydroData(struct hydro_dataframe *hydro_data, int frame, double inj_radius, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr)
{
    //wrapper function that collects the hydro data based on the defined parameters from the user
    char hydro_file[STR_BUFFER]="";
    char hydro_prefix[STR_BUFFER]="";
    
    snprintf(hydro_prefix,sizeof(hydro_prefix),"%s%s",FILEPATH,FILEROOT );

    #if DIMENSIONS == TWO
    
        #if SIM_SWITCH == FLASH
            //if using FLASH data for 2D
            //put proper number at the end of the flash file
            modifyFlashName(hydro_file, hydro_prefix, frame);
            
            fprintf(fPtr,">> MCRaT is opening FLASH file %s\n", hydro_file);
            fflush(fPtr);
            
            readAndDecimate(hydro_file, hydro_data, inj_radius, ph_inj_switch, min_r, max_r, min_theta, max_theta, fPtr);
        #elif SIM_SWITCH == PLUTO_CHOMBO
            modifyPlutoName(hydro_file, hydro_prefix, frame);
            
            fprintf(fPtr,">> MCRaT is opening PLUTO-Chombo file %s\n", hydro_file);
            fflush(fPtr);
            
            readPlutoChombo(hydro_file, hydro_data, inj_radius, ph_inj_switch, min_r, max_r, min_theta, max_theta, fPtr);
        #elif SIM_SWITCH == PLUTO
            modifyPlutoName(hydro_file, hydro_prefix, frame);
            //read in 2D PLUTO data
            fprintf(fPtr,">> MCRaT is opening PLUTO file %s\n", hydro_file);
            fflush(fPtr);
            readPluto(hydro_file, hydro_data, inj_radius, ph_inj_switch, min_r, max_r, min_theta, max_theta, fPtr);
            
        #else
            //if using RIKEN hydro data for 2D szx becomes delta r szy becomes delta theta
            readHydro2D(FILEPATH, frame, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                        &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, ph_inj_switch, min_r, max_r, fPtr);
            //fprintf(fPtr, "%d\n\n", array_num);
        #endif
        
    #else
        #if SIM_SWITCH == FLASH
            #error 3D FLASH simulations are not supported in MCRaT yet.
        #elif SIM_SWITCH == PLUTO_CHOMBO
            modifyPlutoName(hydro_file, hydro_prefix, frame);
            
            fprintf(fPtr,">> MCRaT is opening PLUTO-Chombo file %s\n", hydro_file);
            fflush(fPtr);
            
            readPlutoChombo(hydro_file, hydro_data, inj_radius, ph_inj_switch, min_r, max_r, min_theta, max_theta, fPtr);

        #elif SIM_SWITCH == PLUTO
            modifyPlutoName(hydro_file, hydro_prefix, frame);
            fprintf(fPtr,">> MCRaT is opening PLUTO file %s\n", hydro_file);
            fflush(fPtr);
            readPluto(hydro_file, hydro_data, inj_radius, ph_inj_switch, min_r, max_r, min_theta, max_theta, fPtr);
        #else
            //RKEN Data files
            read_hydro(FILEPATH, frame, inj_radius, &xPtr,  &yPtr, &zPtr,  &szxPtr, &szyPtr, &rPtr,\
                       &thetaPtr, &phiPtr, &velxPtr,  &velyPtr, &velzPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, ph_inj_switch, min_r, max_r, fps_modified, fPtr);
        #endif
    #endif
        
    fprintf(fPtr, "MCRaT: The chosen number of hydro elements is %d\n", hydro_data->num_elements);

    //convert hydro coordinates to spherical so we can inject photons, overwriting values, etc.
    fillHydroCoordinateToSpherical(hydro_data);


    //check for run type see if we need to rewrite any data
    
    #if SIMULATION_TYPE == CYLINDRICAL_OUTFLOW
        cylindricalPrep(hydro_data, fPtr);
    #elif SIMULATION_TYPE == SPHERICAL_OUTFLOW
        sphericalPrep(hydro_data, fPtr);
    #elif SIMULATION_TYPE == STRUCTURED_SPHERICAL_OUTFLOW
        structuredFireballPrep(hydro_data, fPtr);
    #endif

    #if NONTHERMAL_E_DIST != OFF
        calculateElectronDistSubgroupDens(hydro_data->electron_dens_subgroup, fPtr);
        fprintf(fPtr, "electorn dist subgroups: %e %e %e \n", (hydro_data->electron_dens_subgroup)[0], (hydro_data->electron_dens_subgroup)[1], (hydro_data->electron_dens_subgroup)[2]);
        calculateAverageDimlessTheta(hydro_data, fPtr);
        fprintf(fPtr,">> The average dimless temp is %e\n", hydro_data->average_dimless_theta);
        calculateNonthermalElectronDens(hydro_data, fPtr);
    #endif



    return 0;
}

int printHydroGeometry(FILE *fPtr)
{
    
    #if DIMENSIONS == TWO
        char dim[]="2D";
    #elif DIMENSIONS == TWO_POINT_FIVE
        char dim[]="2.5D";
    #else
        char dim[]="3D";
    #endif
    
    #if SIM_SWITCH == FLASH
        char sim[]="Flash";
    #elif SIM_SWITCH == PLUTO_CHOMBO
        char sim[]="PLUTO-Chombo";
    #elif SIM_SWITCH == PLUTO
        char sim[]="PLUTO";
    #endif
    
    #if GEOMETRY == CARTESIAN
        char geo[]="Cartesian";
    #elif GEOMETRY == CYLINDRICAL
        char geo[]="Cylindrical";
    #elif GEOMETRY == SPHERICAL
        char geo[]="Spherical";
    #elif GEOMETRY == POLAR
        char geo[]="Polar";
    #endif
    
    fprintf(fPtr, "MCRaT is working on a %s %s %s simulation.\n", dim, geo, sim );
    fflush(fPtr);
    
    return 0;
}
