/*
  *  This code is to merge all the files among different directories once the MCRaT simulation is complete
  *  The input should be the main CMC directory and the sub directories of each angle range with the number of MPI processes used to start the simulation in each directory
  *  eg call: mpiexec -np X /.merge /dir/to/CMC_dir/  
  *  where X shuould be a multiple of the number of sub directories
  *  SHOULD BE COMPILED WITH -O2 OPTIMIZATION
*/ //TEST WITH PROCESSES INJECTING IN MULTIPLE FRAMES

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "hdf5.h"
#include "mclib.h"
#include "mpi.h"

int main(int argc, char **argv)
{
    double *p0=NULL, *p1=NULL, *p2=NULL, *p3=NULL, *comv_p0=NULL, *comv_p1=NULL, *comv_p2=NULL, *comv_p3=NULL, *r0=NULL, *r1=NULL, *r2=NULL, *s0=NULL, *s1=NULL, *s2=NULL, *s3=NULL, *num_scatt=NULL, *weight=NULL;
    double *p0_p=NULL, *p1_p=NULL, *p2_p=NULL, *p3_p=NULL, *comv_p0_p=NULL, *comv_p1_p=NULL, *comv_p2_p=NULL, *comv_p3_p=NULL, *r0_p=NULL, *r1_p=NULL, *r2_p=NULL, *s0_p=NULL, *s1_p=NULL, *s2_p=NULL, *s3_p=NULL, *num_scatt_p=NULL, *weight_p=NULL;
    int num_angle_dirs=0, i=0, j=0, k=0, l=0, num_types=12;
    int *num_procs_per_dir=NULL, frm0_small, frm0_large, last_frm, frm2_small, frm2_large, small_frm, large_frm, frm=0, all_photons;
    int *frm_array=NULL, *each_subdir_number=NULL, *displPtr=NULL;
    int myid, numprocs, subdir_procs, subdir_id, frames_to_merge, start_count, end_count ;
    int  count=0, index=0,  isNotCorrupted=0;
    int file_count = 0, max_num_procs_per_dir=0;
    int *photon_injection_count=NULL;
    double garbage;
    char mc_file[500]="" ;
    char dir[500]="";
    char group[500]="";
    char merged_filename[500]="";
    char filename_k[2000]="", mcdata_type[20]="";
    char *str="mc_proc_", *ph_type=NULL, *ph_type_p=NULL;
    struct dirent* dent;
    DIR * dirp;
    struct dirent * entry;
    DIR* srcdir = opendir(argv[1]);
    
    MPI_Datatype stype;
    hid_t       file, file_id, group_id, dspace, fspace, mspace;         /* file and dataset identifiers */
    hid_t	plist_id_file, plist_id_data;        /* property list identifier( access template) */
    hsize_t dims[1]={0},dims_old[1]={0};
     hsize_t maxdims[1]={H5S_UNLIMITED};
     hsize_t      size[1];
    hsize_t      offset[1];
    herr_t	status, status_group;
    hid_t dset_p0, dset_p1, dset_p2, dset_p3, dset_comv_p0, dset_comv_p1, dset_comv_p2, dset_comv_p3, dset_r0, dset_r1, dset_r2, dset_s0, dset_s1, dset_s2, dset_s3, dset_num_scatt, dset_weight, dset_ph_type;
    
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


    while((dent = readdir(srcdir)) != NULL)
    {
        struct stat st;
        //file_count=0;
        
        if(strcmp(dent->d_name, ".") == 0 || strcmp(dent->d_name, "..") == 0)
            continue;

        if (fstatat(dirfd(srcdir), dent->d_name, &st, 0) < 0)
        {
            perror(dent->d_name);
            continue;
        }

        if (S_ISDIR(st.st_mode)) 
        {
            //snprintf(dir,sizeof(dir),"%s",dent->d_name );
            if (strstr(dent->d_name, "ALL_DATA") == NULL)
            {
                num_angle_dirs++;
                //printf("found directory %s\n", dent->d_name);
            }
        }
        
    }
    
    closedir(srcdir);
    
    num_procs_per_dir=malloc(num_angle_dirs*sizeof(int));
    each_subdir_number=malloc(num_angle_dirs*sizeof(int));
    displPtr=malloc(num_angle_dirs*sizeof(int));
    *(displPtr+0)=0;
     char *dirs[num_angle_dirs];
    
    count=0;
    srcdir = opendir(argv[1]);
    while((dent = readdir(srcdir)) != NULL)
    {
        struct stat st;
        file_count=0;
        
        if(strcmp(dent->d_name, ".") == 0 || strcmp(dent->d_name, "..") == 0)
            continue;

        if (fstatat(dirfd(srcdir), dent->d_name, &st, 0) < 0)
        {
            perror(dent->d_name);
            continue;
        }

        if (S_ISDIR(st.st_mode)) 
        {
            //printf("found directory %s\n", dent->d_name);
             if (strstr(dent->d_name, "ALL_DATA") == NULL)
            {
                snprintf(dir,sizeof(dir),"%s%s/",argv[1],dent->d_name );
                dirs[count] =  malloc((strlen(dir)+1));
                strcpy(dirs[count],dir);
                //printf("SECOND: found directory %s\n", dirs[count]);
            
                dirp = opendir(dir); //do into the directory to get each file
                while ((entry = readdir(dirp)) != NULL) 
                {
                    if ((entry->d_type == DT_REG) && (strstr(entry->d_name, str) != NULL))
                    { /* If the entry is a regular file  */
                        file_count++;
                        //printf("%s\n", entry->d_name );
                    }
                }
                *(num_procs_per_dir +count)=file_count;
                if (max_num_procs_per_dir<file_count)
                {
                    max_num_procs_per_dir=file_count; //find the max number of processes in each directory
                }
                
                count++;
            }
        }
        
    }
    
    closedir(srcdir);
    
    //find number of directories for each angle range
    //printf("%s: %d\n", argv[1], num_angle_dirs);
    
    //for (i=0;i<num_angle_dirs;i++)
    //{
     //   printf("%d\n", *(num_procs_per_dir+i));
     //   printf(" %s\n",  dirs[i]);
    //}
    
    //get the last and initial hydro file in sim
    snprintf(mc_file,sizeof(mc_file),"%s%s",argv[1],MCPAR);
    readMcPar(mc_file, &garbage, &garbage, &garbage,&garbage, &garbage, &garbage, &garbage,&garbage, &frm0_small,&frm0_large, &last_frm ,&frm2_small, &frm2_large, &garbage, &garbage, &i, &i, &i,&i); //thetas that comes out is in degrees
        //printf("%s frm_0small: %d frm_0large: %d, last: %d\n", mc_file, frm0_small,frm0_large, last_frm);
    
    //with all the info make array of all the files that need to be created
    small_frm= (frm0_small < frm0_large) ? frm0_small : frm0_large;
    large_frm= (frm2_small > frm2_large) ? frm2_small : frm2_large;
    frm_array=malloc(sizeof(int)*(last_frm-small_frm+1));
    count=0;
    for (i=small_frm;i<last_frm+1;i++)
    {
        //printf("Count: %d\n",count);
        *(frm_array+count)=i;
        count++;
    }
    
    //set up the ALL_DATA directory name
    snprintf(dir,sizeof(dir),"%sALL_DATA/",argv[1] );
    
    //set up MPI and break up the processes into groups of the number of sub directories
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    //break up into groups by subdir
    index= myid/num_angle_dirs;
    MPI_Comm frames_to_merge_comm;
    MPI_Comm_split(MPI_COMM_WORLD, index , myid, &frames_to_merge_comm);
    MPI_Comm_rank(frames_to_merge_comm, &subdir_id);
    MPI_Comm_size(frames_to_merge_comm, &subdir_procs);
    
    frames_to_merge=(count-1)/(numprocs/num_angle_dirs); //count-1 b/c @ end of loop it added 1
    start_count=*(frm_array+(index*frames_to_merge));
    end_count=*(frm_array+((index+1)*frames_to_merge));
    

    if (index==(numprocs/num_angle_dirs)-1)
    {
        //printf("in If\n");
        end_count=*(frm_array+(count-1))+1;
    }
    
    
    
    
    printf("subdir_id %d, subdir_procs %d subdir %s, num of frames to merge %d, index %d\n", subdir_id, subdir_procs, dirs[subdir_id], frames_to_merge, index );
    printf("Start file %d end file %d\n", start_count, end_count);
    
    //exit(0);
    
    if (myid==0)
    {
        //have 1st process see if the folder ALL_DATA exists and if not create it
        dirp = opendir(dir);
        if (ENOENT == errno)
        {
            //if it doesnt exist create it
            mkdir(dir, 0777); //make the directory with full permissions
        }
        else
        {
            closedir(dirp);
        }
    }
    
    //directory exists now, can create files in it with appropriate datasets, all processes in communicator participate in this
    //Set up file access property list with parallel I/O access
    MPI_Info info  = MPI_INFO_NULL;
    
    photon_injection_count=malloc((*(num_procs_per_dir+subdir_id))*sizeof(int)); //to incrememnt the number of photons already injected by a process
    for (k=0;k<*(num_procs_per_dir+subdir_id);k++)
    {
        *(photon_injection_count+k)=0;
    }
    
    //create files
    //start_count=2474;
    //end_count=143;
    for (i= end_count-1; i>=start_count;i--)
    {
        //go through the mpi files to find the total number of photons needed for the final dataset  
        //printf("\n\n%d\n", i);
        dims[0]=0;
         j=0;
        for (k=0;k<*(num_procs_per_dir+subdir_id);k++)
        {
            //for each process' file, find out how many elements and add up to find total number of elements needed in the data set for the frame number
            snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dirs[subdir_id],"mc_proc_", k, ".h5" );
            //printf("Dir: %s\n",filename_k );
            
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
                
                //printf("File %s num_ph %d\n", filename_k, j);
                
                status = H5Sclose (dspace);
                status = H5Dclose (dset_p0);
                status = H5Gclose(group_id);
            }
            status = H5Fclose(file);
            
        }
        
        //find total number of photons
        MPI_Allreduce(&j, &all_photons, 1, MPI_INT, MPI_SUM, frames_to_merge_comm);
        
        //get the number for each subdir for later use
        //MPI_Allgather(&j, 1, MPI_INT, each_subdir_number, 1, MPI_INT, frames_to_merge_comm);
        
        //set up the displacement of data
        //for (j=1;j<num_angle_dirs;j++)
        //{
        //    *(displPtr+j)=(*(displPtr+j-1))+(*(each_subdir_number+j-1));
        //}
        
        //if (subdir_id==0)
        //{
        //    printf("Frame: %d Total photons %d\n", i, all_photons);
        //}
        
        
        plist_id_file = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id_file, frames_to_merge_comm, info);
        
        snprintf(merged_filename,sizeof(merged_filename),"%smcdata_%d.h5",dir, i );
        status = H5Eset_auto(NULL, NULL, NULL); //turn off automatic error printing
        file_id = H5Fcreate(merged_filename, H5F_ACC_EXCL, H5P_DEFAULT, plist_id_file);
        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr); //turn on auto error printing
        
        //if the file exists we have to check it to ensure its not corrupted
        
         if (file_id<0)
        {
            //printf( "Checking File %s\n",merged_filename );
            //the file exists, open it with read write 
            file_id=H5Fopen(merged_filename, H5F_ACC_RDWR, plist_id_file);
            
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
                dset_p0 = H5Dopen (file_id, mcdata_type, H5P_DEFAULT); //open dataset
                
                //get the number of points
                dspace = H5Dget_space (dset_p0);
                status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                
                //fprintf(fPtr, "j:%d, dim: %d\n",j, dims[0] );
                //fflush(fPtr);
                
                isNotCorrupted += fmod(dims[0], all_photons); //if the dimension is the dame then the fmod ==0 (remainder of 0), if all datatsets are ==0 then you get a truth value of 0 meaning that it isnt corrupted
                
                status = H5Sclose (dspace);
                status = H5Dclose (dset_p0);
            }
            
            status = H5Fclose(file_id);
            file_id=-1; //do this so if the file exists it doesnt go into the rewriting portion just based on that
        }
        
        //printf("file %s has isNotCorrupted=%d\n", merged_filename, isNotCorrupted );
        
        
        if ((file_id>=0) || (isNotCorrupted != 0 ))
        {
            if (isNotCorrupted != 0)
            {
                //if the data is corrupted overwrite the file
                file_id = H5Fcreate(merged_filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id_file);
            }
            
            frm=(i <  large_frm) ? i  : large_frm;
            
            for (k=0;k<*(num_procs_per_dir+subdir_id);k++)
            {
                *(photon_injection_count+k)=0; //reset the count to 0
            }
            
            //order data based on which processes injected photons 1st
            #if SYNCHROTRON_SWITCH == ON
            l=i;
            #else
            for (l=small_frm;l<frm+1;l++)
            #endif
            {
                //printf("\n\n %d\n\n",l);
                //read the data in from each process in a given subdir, use max_num_procs_per_dir in case one directory used more processes than the others and deal with it in code
                for (k=0;k<max_num_procs_per_dir;k++)
                {                
                
                    dims[0]=0;
                    j=0;
                    //for each process' file, find out how many elements and add up to find total number of elements needed in the data set for the frame number
                    snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dirs[subdir_id],"mc_proc_", k, ".h5" );
                    //printf("Dir: %s\n",filename_k );
            
                    if (k<*(num_procs_per_dir+subdir_id))
                    {
                        //we know that the process exists and the file should exist
                        //open the file
                        status = H5Eset_auto(NULL, NULL, NULL); //turn of error printing if the file doesnt exist, if the process number doesnt exist
                        file=H5Fopen(filename_k, H5F_ACC_RDONLY, H5P_DEFAULT);
                        status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
                    }
                    else
                    {
                        //know that the process doesnt exist within that subdirectory so dont rry to open a non-existant file
                        file=-1;
                    }
                
                    if (file>=0)
                    {
                        {

                            //see if the frame exists
                            /*
                            snprintf(group,sizeof(group),"%d",i );
                            status = H5Eset_auto(NULL, NULL, NULL);
                            status_group = H5Gget_objinfo (file, group, 0, NULL);
                            status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
                            */
                            snprintf(group,sizeof(group),"%d/PW",l );
                            status = H5Eset_auto(NULL, NULL, NULL);
                            status_group = H5Gget_objinfo (file, group, 0, NULL);
                            status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
                        }

                    }
                
            
                    //if it does open it and read in the size
                    //#if SYNCHROTRON_SWITCH == ON
                    //if (status_group >= 0 && file>=0 && l>=i)
                    //#else
                    if (status_group >= 0 && file>=0)
                    //#endif
                    {
                        //read in the number of injected photons first
                        #if SYNCHROTRON_SWITCH == ON
                            snprintf(group,sizeof(group),"%d",i );
                        #else
                            snprintf(group,sizeof(group),"%d",l );
                        #endif
                        group_id = H5Gopen2(file, group, H5P_DEFAULT);
                        dset_weight = H5Dopen (group_id, "PW", H5P_DEFAULT);
                        dspace = H5Dget_space (dset_weight);
                        status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                        j=dims[0];//calculate the total number of photons to save to new hdf5 file
                        status = H5Sclose (dspace);
                        status = H5Dclose (dset_weight);
                        status = H5Gclose(group_id);
                        //printf("Num of ph: %d\n", j);
                        
                        snprintf(group,sizeof(group),"%d",i ); 
                    
                        //printf("Opening dataset\n");
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
                        
                        #if SYNCHROTRON_SWITCH == ON
                        {
                            dset_weight = H5Dopen (group_id, "PW", H5P_DEFAULT);
                        }
                        #else
                        {
                            dset_weight = H5Dopen (file, "PW", H5P_DEFAULT);//for non synch runs look at the global /PW dataset
                        }
                        #endif
                        
                        #if SAVE_TYPE == ON
                        {
                            dset_ph_type = H5Dopen (group_id, "PT", H5P_DEFAULT);
                        }
                        #endif

                        
                        //malloc memory
                        p0_p=malloc(j*sizeof(double));  p1_p=malloc(j*sizeof(double));  p2_p=malloc(j*sizeof(double));  p3_p=malloc(j*sizeof(double));
                        
                        #if COMV_SWITCH == ON
                        {
                            comv_p0_p=malloc(j*sizeof(double));  comv_p1_p=malloc(j*sizeof(double));  comv_p2_p=malloc(j*sizeof(double));  comv_p3_p=malloc(j*sizeof(double));
                        }
                        #endif
                        
                        r0_p=malloc(j*sizeof(double));  r1_p=malloc(j*sizeof(double));  r2_p=malloc(j*sizeof(double));
                        
                        #if STOKES_SWITCH == ON
                        {
                            s0_p=malloc(j*sizeof(double));  s1_p=malloc(j*sizeof(double));  s2_p=malloc(j*sizeof(double));  s3_p=malloc(j*sizeof(double));
                        }
                        #endif
                        
                        #if SAVE_TYPE == ON
                        {
                            ph_type_p=malloc((j)*sizeof(char));
                        }
                        #endif
                        
                        num_scatt_p=malloc(j*sizeof(double));
                        
                        weight_p=malloc(j*sizeof(double));
                        
                        //printf("file %d frame: %d, process  %d start: %d, j: %d\n", i, l, k, *(photon_injection_count+k), dims[0]);
                        
                        #if SYNCHROTRON_SWITCH == ON
                        {
                            offset[0]=0;
                        }
                        #else
                        {
                            offset[0]=*(photon_injection_count+k);
                        }
                        #endif
                        
                        //have to read in the data from *(photon_injection_count+k) to *(photon_injection_count+k)+j
                        mspace = H5Screate_simple (1, dims, NULL);
                        dspace = H5Dget_space(dset_p0);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        
                
                        //read the data in
                        status = H5Dread(dset_p0, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (p0_p));
                        status = H5Sclose (dspace);  status = H5Dclose (dset_p0);
                        
                        dspace = H5Dget_space(dset_p1);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_p1, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (p1_p));
                        status = H5Sclose (dspace);  status = H5Dclose (dset_p1);
                        
                        dspace = H5Dget_space(dset_p2);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_p2, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (p2_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_p2); 
                        
                        dspace = H5Dget_space(dset_p3);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_p3, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (p3_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_p3);
                        
                        #if COMV_SWITCH == ON
                        {
                            dspace = H5Dget_space(dset_comv_p0);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_comv_p0, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (comv_p0_p));
                            status = H5Sclose (dspace);  status = H5Dclose (dset_comv_p0);
                            
                            dspace = H5Dget_space(dset_comv_p1);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_comv_p1, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (comv_p1_p));
                            status = H5Sclose (dspace);  status = H5Dclose (dset_comv_p1);
                            
                            dspace = H5Dget_space(dset_comv_p2);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_comv_p2, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (comv_p2_p));
                            status = H5Sclose (dspace);  status = H5Dclose (dset_comv_p2);
                            
                            dspace = H5Dget_space(dset_comv_p3);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_comv_p3, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (comv_p3_p));
                            status = H5Sclose (dspace);  status = H5Dclose (dset_comv_p3);
                        }
                        #endif
                        
                        dspace = H5Dget_space(dset_r0);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_r0, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (r0_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_r0); 
                        
                        dspace = H5Dget_space(dset_r1);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_r1, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (r1_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_r1); 
                        
                        dspace = H5Dget_space(dset_r2);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_r2, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (r2_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_r2);
                        
                        #if STOKES_SWITCH == ON
                        {
                            dspace = H5Dget_space(dset_s0);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_s0, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (s0_p));
                            status = H5Sclose (dspace); status = H5Dclose (dset_s0);
                            
                            dspace = H5Dget_space(dset_s1);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_s1, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (s1_p));
                            status = H5Sclose (dspace); status = H5Dclose (dset_s1);
                            
                            dspace = H5Dget_space(dset_s2);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_s2, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (s2_p));
                            status = H5Sclose (dspace); status = H5Dclose (dset_s2);
                            
                            dspace = H5Dget_space(dset_s3);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_s3, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (s3_p));
                            status = H5Sclose (dspace); status = H5Dclose (dset_s3);
                        }
                        #endif
                        
                        dspace = H5Dget_space(dset_num_scatt);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_num_scatt, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (num_scatt_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_num_scatt);
                        
                        //printf("Before Weight read\n");
                        dspace = H5Dget_space(dset_weight);
                        status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dread(dset_weight, H5T_NATIVE_DOUBLE, mspace, dspace, H5P_DEFAULT, (weight_p));
                        status = H5Sclose (dspace); status = H5Dclose (dset_weight);
                        //printf("After Weight read\n");
                        
                        #if SAVE_TYPE == ON
                        {
                            dspace = H5Dget_space(dset_ph_type);
                            status = H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dread(dset_ph_type, H5T_NATIVE_CHAR, mspace, dspace, H5P_DEFAULT, (ph_type_p));
                            status = H5Sclose (dspace); status = H5Dclose (dset_ph_type);
                        }
                        #endif
                        
                        status = H5Sclose (mspace);
                        status = H5Gclose(group_id);
                        
                        
                        
                        //#if SYNCHROTRON_SWITCH == ON
                        //{
                        //    *(photon_injection_count+k)+=0;
                        //}
                        //#else
                        {
                            *(photon_injection_count+k)+=j;
                        }
                        //#endif
                    }
                    else
                    {
                        //allocate memory so Allgather doesn't fail with NULL pointer
                        j=1;
                        p0_p=malloc(j*sizeof(double));  p1_p=malloc(j*sizeof(double));  p2_p=malloc(j*sizeof(double));  p3_p=malloc(j*sizeof(double));
                        
                        #if COMV_SWITCH == ON
                        {
                            comv_p0_p=malloc(j*sizeof(double));  comv_p1_p=malloc(j*sizeof(double));  comv_p2_p=malloc(j*sizeof(double));  comv_p3_p=malloc(j*sizeof(double));
                        }
                        #endif
                        
                        r0_p=malloc(j*sizeof(double));  r1_p=malloc(j*sizeof(double));  r2_p=malloc(j*sizeof(double));
                        
                        #if STOKES_SWITCH == ON
                        {
                            s0_p=malloc(j*sizeof(double));  s1_p=malloc(j*sizeof(double));  s2_p=malloc(j*sizeof(double));  s3_p=malloc(j*sizeof(double));
                        }
                        #endif
                        
                        #if SAVE_TYPE == ON
                        {
                            ph_type_p=malloc((j)*sizeof(char));
                        }
                        #endif
                        
                        num_scatt_p=malloc(j*sizeof(double));
                        
                        weight_p=malloc(j*sizeof(double));
                    }
                    
                    //find total number of photons
                    MPI_Allreduce(&dims[0], &all_photons, 1, MPI_INT, MPI_SUM,  frames_to_merge_comm);
                    //dims[0]=all_photons;
                    
                    //printf("ID %d j: %d\n", subdir_id, dims[0]);
                    
                    //get the number for each subdir for later use
                    MPI_Allgather(&dims[0], 1, MPI_INT, each_subdir_number, 1, MPI_INT,   frames_to_merge_comm);
                    //for (j=0;j<num_angle_dirs;j++)
                    //{
                     //    printf("ID %d eachsubdir_num %d \n",  subdir_id, *(each_subdir_number+j));
                    //}
                    
        
                    //set up the displacement of data
                    for (j=1;j<num_angle_dirs;j++)
                    {
                        *(displPtr+j)=(*(displPtr+j-1))+(*(each_subdir_number+j-1));
                         //printf("Displ %d eachsubdir_num %d \n",  *(displPtr+j), *(each_subdir_number+j-1));
                    }
        
                    //if (subdir_id==0)
                    //{
                    //    printf("Frame: %d Total photons %d\n", i, all_photons);
                    //}
                    
                    //now allocate enough ememory for all_photons in the mpi files from proc 0 initially 
                    p0=malloc(all_photons*sizeof(double));  p1=malloc(all_photons*sizeof(double));  p2=malloc(all_photons*sizeof(double));  p3=malloc(all_photons*sizeof(double));
                    
                    #if COMV_SWITCH == ON
                    {
                        comv_p0=malloc(all_photons*sizeof(double));  comv_p1=malloc(all_photons*sizeof(double));  comv_p2=malloc(all_photons*sizeof(double));  comv_p3=malloc(all_photons*sizeof(double));
                    }
                    #endif
                    
                    r0=malloc(all_photons*sizeof(double));  r1=malloc(all_photons*sizeof(double));  r2=malloc(all_photons*sizeof(double));
                    
                    #if STOKES_SWITCH == ON
                    {
                        s0=malloc(all_photons*sizeof(double));  s1=malloc(all_photons*sizeof(double));  s2=malloc(all_photons*sizeof(double));  s3=malloc(all_photons*sizeof(double));
                    }
                    #endif
                    
                    #if SAVE_TYPE == ON
                    {
                        ph_type=malloc(all_photons*sizeof(char));
                    }
                    #endif
                    
                    num_scatt=malloc(all_photons*sizeof(double));
                    
                    weight=malloc(all_photons*sizeof(double));
                    
                    
                    //save data in correct order to p0, s0, r0, etc. in order of angle 
                    
                    //MPI_Type_commit( &stype ); 
                    MPI_Allgatherv(p0_p, dims[0], MPI_DOUBLE, p0, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    MPI_Allgatherv(p1_p, dims[0], MPI_DOUBLE, p1, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    MPI_Allgatherv(p2_p, dims[0], MPI_DOUBLE, p2, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    MPI_Allgatherv(p3_p, dims[0], MPI_DOUBLE, p3, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    
                    #if COMV_SWITCH == ON
                    {
                        MPI_Allgatherv(comv_p0_p, dims[0], MPI_DOUBLE, comv_p0, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                        MPI_Allgatherv(comv_p1_p, dims[0], MPI_DOUBLE, comv_p1, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                        MPI_Allgatherv(comv_p2_p, dims[0], MPI_DOUBLE, comv_p2, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                        MPI_Allgatherv(comv_p3_p, dims[0], MPI_DOUBLE, comv_p3, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    }
                    #endif
                    
                    MPI_Allgatherv(r0_p, dims[0], MPI_DOUBLE, r0, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    MPI_Allgatherv(r1_p, dims[0], MPI_DOUBLE, r1, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    MPI_Allgatherv(r2_p, dims[0], MPI_DOUBLE, r2, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    
                    #if STOKES_SWITCH == ON
                    {
                        MPI_Allgatherv(s0_p, dims[0], MPI_DOUBLE, s0, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                        MPI_Allgatherv(s1_p, dims[0], MPI_DOUBLE, s1, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                        MPI_Allgatherv(s2_p, dims[0], MPI_DOUBLE, s2, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                        MPI_Allgatherv(s3_p, dims[0], MPI_DOUBLE, s3, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    }
                    #endif
                    
                    #if SAVE_TYPE == ON
                    {
                        MPI_Allgatherv(ph_type_p, dims[0], MPI_CHAR, ph_type, each_subdir_number, displPtr, MPI_CHAR, frames_to_merge_comm);
                    }
                    #endif

                    
                    MPI_Allgatherv(num_scatt_p, dims[0], MPI_DOUBLE, num_scatt, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    
                    MPI_Allgatherv(weight_p, dims[0], MPI_DOUBLE, weight, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                    
                    /*
                    if (subdir_id==0)
                    {
                        for (j=0;j<all_photons;j++)
                        {
                            printf("Read Data: %e Gathered data: %e\n", *(s0+j), *(s0_p+j));
                        }
                    }
                     */
                    //exit(0);
                    
                    dims[0]=all_photons;
                    //if ((k==0)  && (all_photons>0))
                    #if SYNCHROTRON_SWITCH == ON
                    if ((l==i) && (k==0) && (all_photons>0))
                    #else
                    if ((l==small_frm)  && (all_photons>0))
                    #endif
                    {
                        //printf("IN THE IF STATEMENT\n");
                        //set up new dataset
                        //create the datasets with the appropriate number of elements
                        
                        plist_id_data = H5Pcreate (H5P_DATASET_CREATE);
                        status = H5Pset_chunk (plist_id_data, 1, dims);
                        dspace = H5Screate_simple (1, dims, maxdims);
                        
                        
                        
                        dset_p0=H5Dcreate2(file_id, "P0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        dset_p1=H5Dcreate2(file_id, "P1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        dset_p2=H5Dcreate2(file_id, "P2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        dset_p3=H5Dcreate2(file_id, "P3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        
                        #if COMV_SWITCH == ON
                        {
                            dset_comv_p0=H5Dcreate2(file_id, "COMV_P0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                            dset_comv_p1=H5Dcreate2(file_id, "COMV_P1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                            dset_comv_p2=H5Dcreate2(file_id, "COMV_P2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                            dset_comv_p3=H5Dcreate2(file_id, "COMV_P3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        }
                        #endif
                        
                        dset_r0=H5Dcreate2(file_id, "R0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        dset_r1=H5Dcreate2(file_id, "R1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        dset_r2=H5Dcreate2(file_id, "R2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        
                        #if STOKES_SWITCH == ON
                        {
                            dset_s0=H5Dcreate2(file_id, "S0", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                            dset_s1=H5Dcreate2(file_id, "S1", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                            dset_s2=H5Dcreate2(file_id, "S2", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                            dset_s3=H5Dcreate2(file_id, "S3", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        }
                        #endif
                        
                        #if SAVE_TYPE == ON
                        {
                            dset_ph_type=H5Dcreate2(file_id, "PT", H5T_NATIVE_CHAR, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        }
                        #endif

                        
                        dset_num_scatt=H5Dcreate2(file_id, "NS", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                        
                        dset_weight=H5Dcreate2(file_id, "PW", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);

                                                 
                        H5Pclose(plist_id_data);
                        H5Sclose(dspace);
                        
                        plist_id_data = H5Pcreate (H5P_DATASET_XFER);
                        H5Pset_dxpl_mpio (plist_id_data, H5FD_MPIO_COLLECTIVE);
                        
                        //write data
                        offset[0]=0;
                        dspace = H5Dget_space(dset_p0);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, p0);
                        H5Sclose(dspace);
                        
                        
                         dspace = H5Dget_space(dset_p1);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, p1);
                        H5Sclose(dspace);
                        
                        dspace = H5Dget_space(dset_p2);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, p2);
                        H5Sclose(dspace);
                        
                        dspace = H5Dget_space(dset_p3);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, p3);
                        H5Sclose(dspace);
                        
                        #if COMV_SWITCH == ON
                        {
                            dspace = H5Dget_space(dset_comv_p0);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_comv_p0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, comv_p0);
                            H5Sclose(dspace);
                            
                            
                            dspace = H5Dget_space(dset_comv_p1);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_comv_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, comv_p1);
                            H5Sclose(dspace);
                            
                            dspace = H5Dget_space(dset_comv_p2);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_comv_p2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, comv_p2);
                            H5Sclose(dspace);
                            
                            dspace = H5Dget_space(dset_comv_p3);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_comv_p3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, comv_p3);
                            H5Sclose(dspace);
                        }
                        #endif
                        
                        dspace = H5Dget_space(dset_r0);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_r0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, r0);
                        H5Sclose(dspace);
        
                        dspace = H5Dget_space(dset_r1);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_r1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, r1);
                        H5Sclose(dspace);
                        
                        dspace = H5Dget_space(dset_r2);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_r2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, r2);
                        H5Sclose(dspace);
                        
                        #if STOKES_SWITCH == ON
                        {
                            dspace = H5Dget_space(dset_s0);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_s0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, s0);
                            H5Sclose(dspace);
                            
                            dspace = H5Dget_space(dset_s1);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_s1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, s1);
                            H5Sclose(dspace);
                            
                            dspace = H5Dget_space(dset_s2);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_s2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, s2);
                            H5Sclose(dspace);
                            
                            dspace = H5Dget_space(dset_s3);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_s3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  plist_id_data, s3);
                            H5Sclose(dspace);
                        }
                        #endif
                        
                        #if SAVE_TYPE == ON
                        {
                            dspace = H5Dget_space(dset_ph_type);
                            status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            status = H5Dwrite (dset_ph_type, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id_data, ph_type);
                            H5Sclose(dspace);

                        }
                        #endif

                        
                        dspace = H5Dget_space(dset_num_scatt);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, num_scatt);
                        H5Sclose(dspace);
                        
                        dspace = H5Dget_space(dset_weight);
                        status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, weight);
                        H5Sclose(dspace);

                        
                        status = H5Dclose (dset_p0); 
                        status = H5Dclose (dset_p1); status = H5Dclose (dset_p2); status = H5Dclose (dset_p3);
                        
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
                        
                        H5Pclose(plist_id_data);
                        
                    }
                    #if SYNCHROTRON_SWITCH == ON
                    else
                    #else
                    else if ((l>small_frm) && (all_photons>0))
                    #endif
                    {
                        //printf("IN THE ELSE IF STATEMENT\n"); if ((k>0)  && (all_photons>0))
                        plist_id_data = H5Pcreate (H5P_DATASET_XFER);
                        H5Pset_dxpl_mpio (plist_id_data, H5FD_MPIO_COLLECTIVE);
                        
                        dset_p0 = H5Dopen (file_id, "P0", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_p0);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_p0, size);
                        
                        fspace = H5Dget_space (dset_p0);
                        offset[0] = dims_old[0];
                        
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_p0, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, p0);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_p0);
                        
                        dset_p1 = H5Dopen (file_id, "P1", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_p1);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_p1, size);
                        fspace = H5Dget_space (dset_p1);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_p1, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, p1);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_p1);
                        
                        dset_p2 = H5Dopen (file_id, "P2", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_p2);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_p2, size);
                        fspace = H5Dget_space (dset_p2);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_p2, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, p2);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_p2);
                        
                        dset_p3 = H5Dopen (file_id, "P3", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_p3);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_p3, size);
                        fspace = H5Dget_space (dset_p3);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_p3, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, p3);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_p3);
                        
                        #if COMV_SWITCH == ON
                        {
                            dset_comv_p0 = H5Dopen (file_id, "COMV_P0", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_comv_p0);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_comv_p0, size);
                            fspace = H5Dget_space (dset_comv_p0);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_comv_p0, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, comv_p0);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_comv_p0);
                            
                            dset_comv_p1 = H5Dopen (file_id, "COMV_P1", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_comv_p1);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_comv_p1, size);
                            fspace = H5Dget_space (dset_comv_p1);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_comv_p1, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, comv_p1);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_comv_p1);
                            
                            dset_comv_p2 = H5Dopen (file_id, "COMV_P2", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_comv_p2);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_comv_p2, size);
                            fspace = H5Dget_space (dset_comv_p2);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_comv_p2, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, comv_p2);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_comv_p2);
                            
                            dset_comv_p3 = H5Dopen (file_id, "COMV_P3", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_comv_p3);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_comv_p3, size);
                            fspace = H5Dget_space (dset_comv_p3);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_comv_p3, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, comv_p3);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_comv_p3);
                        }
                        #endif
                        
                        dset_r0 = H5Dopen (file_id, "R0", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_r0);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_r0, size);
                        fspace = H5Dget_space (dset_r0);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_r0, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, r0);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_r0);
                        
                        dset_r1 = H5Dopen (file_id, "R1", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_r1);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_r1, size);
                        fspace = H5Dget_space (dset_r1);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_r1, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, r1);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_r1);
                        
                        dset_r2 = H5Dopen (file_id, "R2", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_r2);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_r2, size);
                        fspace = H5Dget_space (dset_r2);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_r2, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, r2);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_r2);
                        
                        #if STOKES_SWITCH == ON
                        {
                            dset_s0 = H5Dopen (file_id, "S0", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_s0);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_s0, size);
                            fspace = H5Dget_space (dset_s0);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_s0, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, s0);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_s0);
                            
                            dset_s1 = H5Dopen (file_id, "S1", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_s1);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_s1, size);
                            fspace = H5Dget_space (dset_s1);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_s1, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, s1);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_s1);
                            
                            dset_s2 = H5Dopen (file_id, "S2", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_s2);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_s2, size);
                            fspace = H5Dget_space (dset_s2);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_s2, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, s2);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_s2);
                            
                            dset_s3 = H5Dopen (file_id, "S3", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_s3);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_s3, size);
                            fspace = H5Dget_space (dset_s3);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_s3, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, s3);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_s3);
                        }
                        #endif
                        
                        #if SAVE_TYPE == ON
                        {
                            dset_ph_type = H5Dopen (file_id, "PT", H5P_DEFAULT); //open dataset
                            dspace = H5Dget_space (dset_ph_type);
                            status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                            size[0] = dims[0]+ dims_old[0];
                            status = H5Dset_extent (dset_ph_type, size);
                            fspace = H5Dget_space (dset_ph_type);
                            offset[0] = dims_old[0];
                            status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                            mspace = H5Screate_simple (1, dims, NULL);
                            status = H5Dwrite (dset_ph_type, H5T_NATIVE_CHAR, mspace, fspace, plist_id_data, ph_type);
                            status = H5Sclose (dspace);
                            status = H5Sclose (mspace);
                            status = H5Sclose (fspace);
                            status = H5Dclose (dset_ph_type);

                        }
                        #endif

                        
                        dset_num_scatt = H5Dopen (file_id, "NS", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_num_scatt);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_num_scatt, size);
                        fspace = H5Dget_space (dset_num_scatt);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_num_scatt, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, num_scatt);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_num_scatt);
                        
                        //printf("Before weight write\n");
                        dset_weight = H5Dopen (file_id, "PW", H5P_DEFAULT); //open dataset
                        dspace = H5Dget_space (dset_weight);
                        status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        size[0] = dims[0]+ dims_old[0];
                        status = H5Dset_extent (dset_weight, size);
                        fspace = H5Dget_space (dset_weight);
                        offset[0] = dims_old[0];
                        status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                        mspace = H5Screate_simple (1, dims, NULL);
                        status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, weight);
                        status = H5Sclose (dspace);
                        status = H5Sclose (mspace);
                        status = H5Sclose (fspace);
                        status = H5Dclose (dset_weight);
                        //printf("After weight write\n");
                        
                        
                        H5Pclose(plist_id_data);
                        
                    }
                    
                    
                    
                    if (file>=0)
                    {
                        status = H5Fclose(file);
                    }
                    
                    free(p0_p);free(p1_p); free(p2_p);free(p3_p);
                    
                    #if COMV_SWITCH == ON
                    {
                        free(comv_p0_p);free(comv_p1_p); free(comv_p2_p);free(comv_p3_p);
                    }
                    #endif
                    
                    free(r0_p);free(r1_p); free(r2_p);
                    
                    #if STOKES_SWITCH == ON
                    {
                        free(s0_p);free(s1_p); free(s2_p);free(s3_p);
                    }
                    #endif
                    
                    free(num_scatt_p);
                    
                    free(weight_p);
                    
                    #if SAVE_TYPE == ON
                    {
                        free(ph_type_p);
                    }
                    #endif

                    
                    free(p0);free(p1); free(p2);free(p3);
                    
                    #if COMV_SWITCH == ON
                    {
                        free(comv_p0);free(comv_p1); free(comv_p2);free(comv_p3);
                    }
                    #endif
                    
                    free(r0);free(r1); free(r2);
                    
                    #if STOKES_SWITCH == ON
                    {
                        free(s0);free(s1); free(s2);free(s3);
                    }
                    #endif
                    
                    free(num_scatt);
                    
                    free(weight);
                    
                    #if SAVE_TYPE == ON
                    {
                        free(ph_type);
                    }
                    #endif
                    //exit(0);
                }
            
            
            
            }
            H5Fclose(file_id);
        }
        
        
        
        
        H5Pclose(plist_id_file);
        
    }
    
    /*
    if (index==0)
    {
        plist_id_file = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id_file, frames_to_merge_comm, info);
        //merge the weights in the same order
        snprintf(merged_filename,sizeof(merged_filename),"%smcdata_PW.h5",dir );
        
        file_id = H5Fcreate(merged_filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id_file);
        
         for (i= small_frm; i<large_frm+1;i++)
         {
        //last_frm
            for (k=0;k<max_num_procs_per_dir;k++)
            {
                dims[0]=0;
                j=0;
                
                snprintf(filename_k,sizeof(filename_k),"%s%s%d%s",dirs[subdir_id],"mc_proc_", k, ".h5" );
                //printf("Dir: %s\n",filename_k );
            
                //open the file and see if t exists
                status = H5Eset_auto(NULL, NULL, NULL); //turn of error printing if the file doesnt exist, if the process number doesnt exist
                file=H5Fopen(filename_k, H5F_ACC_RDONLY, H5P_DEFAULT);
                status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
                
                if (file>=0)
                {
                    //if the file exists, see if the frame exists
                    snprintf(group,sizeof(group),"%d/PW",i );
                    status = H5Eset_auto(NULL, NULL, NULL);
                    status_group = H5Gget_objinfo (file, group, 0, NULL);
                    status = H5Eset_auto(H5E_DEFAULT, H5Eprint2, stderr);
                }
                
                //printf("Proc %d has status_group %d\n", subdir_id, status_group);
                
                if ((status_group == 0) && (file>=0))
                {
            
                    //read dataset and then 
                     snprintf(group,sizeof(group),"%d",i );
                    group_id = H5Gopen2(file, group, H5P_DEFAULT);
                    dset_weight = H5Dopen (group_id, "PW", H5P_DEFAULT); //open dataset
                    
                    //get the number of points
                    dspace = H5Dget_space (dset_weight);
                    status=H5Sget_simple_extent_dims(dspace, dims, NULL); //save dimesnions in dims
                    j=dims[0];//calculate the total number of photons to save to new hdf5 file
                    
                    weight_p=malloc(j*sizeof(double));
                    
                    status = H5Dread(dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (weight_p));
                    
                    status = H5Sclose (dspace);
                    status = H5Dclose (dset_weight);
                    status = H5Gclose(group_id);
                }
                else
                {
                    //theres nothing to read
                    j=1;
                    weight_p=malloc(j*sizeof(double));
                }
                
                //find total number of photons
                MPI_Allreduce(&dims[0], &all_photons, 1, MPI_INT, MPI_SUM,  frames_to_merge_comm);
                    
                //printf("ID %d j: %d\n", subdir_id, dims[0]);
                    
                //get the number for each subdir for later use
                MPI_Allgather(&dims[0], 1, MPI_INT, each_subdir_number, 1, MPI_INT,   frames_to_merge_comm);
                //for (j=0;j<num_angle_dirs;j++)
                //{
                 //       printf("ID %d eachsubdir_num %d \n",  subdir_id, *(each_subdir_number+j));
                //}
                    
        
                //set up the displacement of data
                for (j=1;j<num_angle_dirs;j++)
                {
                    *(displPtr+j)=(*(displPtr+j-1))+(*(each_subdir_number+j-1));
                       // printf("Displ %d eachsubdir_num %d \n",  *(displPtr+j), *(each_subdir_number+j-1));
                }
                
                weight=malloc(all_photons*sizeof(double)); 
                    
                //MPI_Type_commit( &stype ); 
                MPI_Allgatherv(weight_p, dims[0], MPI_DOUBLE, weight, each_subdir_number, displPtr, MPI_DOUBLE, frames_to_merge_comm);
                
                dims[0]=all_photons;
                
                if ((i== small_frm) && (all_photons>0))
                {
                    //create new datatset
                    plist_id_data = H5Pcreate (H5P_DATASET_CREATE);
                    status = H5Pset_chunk (plist_id_data, 1, dims);
                    dspace = H5Screate_simple (1, dims, maxdims);
                        
                    dset_weight=H5Dcreate(file_id, "PW", H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, plist_id_data, H5P_DEFAULT);
                    H5Pclose(plist_id_data);
                    H5Sclose(dspace);
                        
                    plist_id_data = H5Pcreate (H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio (plist_id_data, H5FD_MPIO_COLLECTIVE);
                        
                    //write data
                    offset[0]=0;
                    dspace = H5Dget_space(dset_weight);
                    status = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                    status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id_data, weight);
                    H5Sclose(dspace);
                    status = H5Dclose (dset_weight);
                        
                    H5Pclose(plist_id_data);
                    
                }
                else if ((i != small_frm) && (all_photons>0))
                {
                    //extend the datatset
                    plist_id_data = H5Pcreate (H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio (plist_id_data, H5FD_MPIO_COLLECTIVE);
                        
                    dset_weight = H5Dopen (file_id, "PW", H5P_DEFAULT); //open dataset
                    dspace = H5Dget_space (dset_weight);
                    status=H5Sget_simple_extent_dims(dspace, dims_old, NULL); //save dimesnions in dims_old
                        
                    size[0] = dims[0]+ dims_old[0];
                    status = H5Dset_extent (dset_weight, size);
                        
                    fspace = H5Dget_space (dset_weight);
                    offset[0] = dims_old[0];
                        
                    status = H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, dims, NULL);
                    mspace = H5Screate_simple (1, dims, NULL);
                    status = H5Dwrite (dset_weight, H5T_NATIVE_DOUBLE, mspace, fspace, plist_id_data, weight);
                    status = H5Sclose (dspace);
                    status = H5Sclose (mspace);
                    status = H5Sclose (fspace);
                    status = H5Dclose (dset_weight);
                }
                
                //exit(0);
                 if (file>=0)
                {
                    status = H5Fclose(file);
                }
            
            }
            
         }
        
        H5Pclose(plist_id_file);
        free(weight_p); free(weight);
        H5Fclose(file_id);
    }
    */
    
    MPI_Finalize();
    
    free(num_procs_per_dir);
    free(each_subdir_number);
    free(displPtr);
    free(frm_array);
}


