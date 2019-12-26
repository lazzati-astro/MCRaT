/*
# Program to run a Monte Carlo radiation transfer through the 2D
# simulations of GRB jets.
#
# Python code written by D. Lazzati at Oregonstate, C code written by Tyler Parsotan @ Oregon State 
# ver 0.1 July 8, 2015
# ver 1.1 July 20, 2015: added record of number of scatterings, included
# 	all terms in weight. Should now give correct light curves.
# ver 1.2 July 21, 2015: added parameter file to keep track of input
# 	params of each simulation

# ver 2.0 July 22, 2015: corrected the problem that arises when there is
# 	no scattering in the time span of one frame. Fixed output arrays dimension.

# ver 2.1 July 25, 2015: fixed bug that did not make the number of
# 	scattering grow with the number of photons.

# ver 3.0 July 28, 2015: using scipy nearest neighbor interpolation to
# 	speed things up. Gained about factor 2

# ver 3.1 July 29, 2015: added radial spread of photon injection points
# ver 3.2 July 31, 2015: added Gamma to the weight of photons!!!

# ver 4.0 Aug 5, 2015: try to speed up by inverting cycle
# ver 4.1 Aug 8, 2015: add spherical test as an option
# ver 4.2 Aug 9, 2015: saving files appending rather than re-writing
# ver 4.3 Aug 11, 2015: corrected error in the calculation of the local temperature
# ver 4.4 Aug 13, 2015: added cylindrical test
# ver 4.5 Aug 18, 2015: fixd various problems pointed by the cylindrical test
# ver 4.6 Aug 21, 2015: corrected mean free path for large radii

# ver 5.0 Aug 25, 2015: corrected problem with high-T electrons and excess scatterings
# ver 5.1 Aug 25, 2015: cleaned-up coding
# ver 5.2 Sept 3, 2015: fixed problem with number of scatterings for multiple injections
 * 
 * ver 6.0 Dec 28, 2016: rewrote the code in C, added checkpoint file so if the code is interrupted all the progress wont be lost, made the code only need to be compiled once for a given MC_XXX directory path
                                            so you just need to supply the sub directory of MC_XXX as a command line argument
*  version 7.0 used OpenMP to parallelize the code by angle and the function findminmfp()
 
 version 8.0 added 3D capabilities for RIKEN hydro data  and 2D capablities for RIKEN 2D hydro data and made it more efficient with grid selection to speed it up
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "mclib.h"
//#include "mclib_3d.h"
//#include "mclib_pluto.h"
#include <omp.h>
#include "mpi.h"

int main(int argc, char **argv)
{
    //compile each time a macro is changed, have to supply the subfolder within the MC_PATH directory as a command line argument to the C program eg. MCRAT 1/
    
	// Define variables
	char flash_prefix[200]="";
	char mc_file[200]="" ;
    char spect;//type of spectrum
    char restrt;//restart or not
    double fps, fps_modified, theta_jmin, theta_jmax,hydro_domain_y, hydro_domain_x ;//frames per second of sim, min opening angle of jet, max opening angle of jet in radians, max y value in hydro domain
    double inj_radius_small, inj_radius_large,  ph_weight_suggest, ph_weight_small, ph_weight_large ;//radius at chich photons are injected into sim
    int frm0,last_frm, frm2_small, frm2_large, j=0, min_photons, max_photons, frm0_small, frm0_large ;//frame starting from, last frame of sim, frame of last injection
    int dim_switch=0;
    int find_nearest_grid_switch=0;
    int increment_inj=1, increment_scatt=1; //increments for injection loop and scattering loop, outer and inner loops respectively, the increment can change for RIKEN 3D hydro files
    
    double inj_radius;
    int frm2,save_chkpt_success=0;
    char mc_filename[200]="";
    char mc_filename_2[200]="";
    char mc_operation[200]="";
    char mc_dir[200]="" ;
    int file_count = 0;
    DIR * dirp;
    struct dirent * entry;
    struct stat st = {0};
    double theta_jmin_thread=0, theta_jmax_thread=0;
        
    char flash_file[200]="";
    char log_file[200]="";
    FILE *fPtr=NULL; //pointer to log file for each thread
    double *xPtr=NULL,  *yPtr=NULL,  *rPtr=NULL,  *thetaPtr=NULL,  *velxPtr=NULL,  *velyPtr=NULL,  *densPtr=NULL,  *presPtr=NULL,  *gammaPtr=NULL,  *dens_labPtr=NULL;
    double *szxPtr=NULL,*szyPtr=NULL, *tempPtr=NULL; //pointers to hold data from FLASH files
    double *phiPtr=NULL, *velzPtr=NULL, *zPtr=NULL, *all_time_steps=NULL;
    int num_ph=0, array_num=0, ph_scatt_index=0, max_scatt=0, min_scatt=0,i=0; //number of photons produced in injection algorithm, number of array elleemnts from reading FLASH file, index of photon whch does scattering, generic counter
    double dt_max=0, thescatt=0, accum_time=0; 
    double  gamma_infinity=0, time_now=0, time_step=0, avg_scatt=0,avg_r=0; //gamma_infinity not used?
    double ph_dens_labPtr=0, ph_vxPtr=0, ph_vyPtr=0, ph_tempPtr=0, ph_vzPtr=0;// *ph_cosanglePtr=NULL ;
    double min_r=0, max_r=0, min_theta=0, max_theta=0;
    int frame=0, scatt_frame=0, frame_scatt_cnt=0, scatt_framestart=0, framestart=0;
    struct photon *phPtr=NULL; //pointer to array of photons 
    
    int angle_count=0;
    int num_angles=0, old_num_angle_procs=0; //old_num_angle_procs is to hold the old number of procs in each angle when cont sims, if  restarting sims this gets set to angle_procs
    int *frame_array=NULL, *proc_frame_array=NULL, *element_num=NULL, *sorted_indexes=NULL, proc_frame_size=0;
    double *thread_theta=NULL; //saves ranges of thetas for each thread to go through
    double delta_theta=1;
    
    int myid, numprocs, angle_procs, angle_id, procs_per_angle;
    
   
   //new OpenMPI stuff
   MPI_Init(NULL,NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   
  
    //new muliple threads injecting and propagating photons
    const gsl_rng_type *rng_t;
    gsl_rng *rng;
    gsl_rng_env_setup();
    rng_t = gsl_rng_ranlxs0;
    rng = gsl_rng_alloc (rng_t); //initalize random number generator to seed the others with random numbers

    
    //want to break up simulation by angle and injection frame & have each thread save data in its own folder 
    //have each thread check if its directory is made and if its restarting (delete evrything) or if its continuing with a previous simulation
    //the angle and the injection frames will be the names of mc_dir, therefore read mc.par first in MC_XXX directory
    
    //make strings of proper directories etc.
	snprintf(flash_prefix,sizeof(flash_prefix),"%s%s",FILEPATH,FILEROOT );
    snprintf(mc_file,sizeof(flash_prefix),"%s%s%s",FILEPATH, MC_PATH,MCPAR);

    
    printf(">> MCRaT:  Reading mc.par: %s\n", mc_file);
    
    readMcPar(mc_file, &hydro_domain_x, &hydro_domain_y, &fps, &theta_jmin, &theta_jmax, &delta_theta, &inj_radius_small,&inj_radius_large, &frm0_small,&frm0_large, &last_frm ,&frm2_small, &frm2_large, &ph_weight_small, &ph_weight_large, &min_photons, &max_photons, &spect, &restrt); //thetas that comes out is in degrees
    //printf("%c\n", restrt);
    
    //divide up angles and frame injections among threads DONT WANT NUMBER OF THREADS TO BE ODD
    //assign ranges to array that hold them
    
    //leave angles in degrees here
    num_angles=(int) (((theta_jmax-theta_jmin)/delta_theta)) ;//*(180/M_PI));
    thread_theta=malloc( num_angles *sizeof(double) );
    *(thread_theta+0)=theta_jmin;//*(180/M_PI);
    //printf("%e\n", *(thread_theta+0));
    
    for (j=1;j<(num_angles); j++)
    {
        *(thread_theta+j)=*(thread_theta+(j-1))+delta_theta;
        //printf("%e\n", *(thread_theta+j));
    }

    

    //make comm without the procs that deal with angle
     //comm for angles
     
     procs_per_angle= numprocs/num_angles;
     //printf("%d\n", procs_per_angle);
     
     MPI_Comm angle_comm;
     if (restrt=='r') //uncomment this when I run MCRAT for sims that didnt originally save angle_procs 
     {
        
        MPI_Comm_split(MPI_COMM_WORLD, myid/procs_per_angle , myid, &angle_comm);
        MPI_Comm_rank(angle_comm, &angle_id);
        MPI_Comm_size(angle_comm, &angle_procs);
        
        //printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myid, numprocs, angle_id, angle_procs);    
                
        theta_jmin_thread= (*(thread_theta+  (myid/procs_per_angle))) *(M_PI/180);
        theta_jmax_thread= theta_jmin_thread+(delta_theta*(M_PI/180));
        
        snprintf(mc_dir,sizeof(flash_prefix),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this
        
        old_num_angle_procs=angle_procs;
     }
     else
     {
         MPI_Group sub_world_group;
         MPI_Comm sub_world_comm;
         int incl_procs[procs_per_angle*num_angles], count, sub_world_id;
         int total_num_to_restart=0;
         int color=1;
         int  *all_cont_process_idPtr=NULL, *each_num_to_restart_per_anglePtr=NULL, *tmp=NULL;
        //for restart='c' case if the number of processes isnt a multiple of procs_per_angle*num_angles make a comm out of those that are in order to analyze files and count number of processes for each angle range need to con't
        count=0;
        for (j=0;j<numprocs;j++)
        {
            if (j<procs_per_angle*num_angles)
            {
                incl_procs[count]=j;
                count++;
            }
        }
        
        if (myid<procs_per_angle*num_angles)
        {
            int myid_2=0;
            // Get the group of processes in MPI_COMM_WORLD and make a sub group to go through checkpoint files
            MPI_Group world_group;
            MPI_Comm root_angle_comm;
            
            MPI_Comm_group(MPI_COMM_WORLD, &world_group);
            MPI_Group_incl(world_group, procs_per_angle*num_angles, incl_procs, &sub_world_group);
            MPI_Comm_create_group(MPI_COMM_WORLD, sub_world_group, 0, &sub_world_comm);
            MPI_Comm_rank(sub_world_comm, &myid_2);
        
            MPI_Comm_split(sub_world_comm, myid_2/procs_per_angle , myid_2, &angle_comm);
            MPI_Comm_rank(angle_comm, &angle_id);
            MPI_Comm_size(angle_comm, &angle_procs);
    
            //create group of all the processes that have angle_id==0
            if (angle_id==0)
            {
                color=0; //set different color for root processes in each group of angle_comm
            }
            MPI_Comm_split(sub_world_comm, color , myid_2, &root_angle_comm); //create comm to exchange info about number of processes to restart for each angle range
            
        
            printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myid, numprocs, angle_id, angle_procs);    
                
            theta_jmin_thread= (*(thread_theta+  (myid_2/procs_per_angle))) *(M_PI/180);
            theta_jmax_thread= theta_jmin_thread+(delta_theta*(M_PI/180));
        
            snprintf(mc_dir,sizeof(flash_prefix),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this
        
            //call the function to count the num of processes for each angle range that need to be con't
            int  count_cont_procs=0, total_cont_procs_angle=0, global_cont_procs=0;
            int *cont_proc_idsPtr=NULL, *total_cont_procs_angle_Ptr=NULL, *displPtr=NULL; //becomes the size of the number of old procceses 
            int *cont_proc_ids_anglePtr=NULL;
            
            old_num_angle_procs=getOrigNumProcesses(&count_cont_procs,  &cont_proc_idsPtr, mc_dir, angle_id,  angle_procs,  last_frm);
            
            if (old_num_angle_procs==-1)
            {
                printf("MCRAT wasnt able to get a value of old_num_angle_procs to continue the simulation. Now exiting to prevent data corruption.\n" );
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            total_cont_procs_angle_Ptr=malloc(angle_procs*sizeof(int));
            displPtr=malloc(angle_procs*sizeof(int));
            MPI_Gather(&count_cont_procs,1,MPI_INT, total_cont_procs_angle_Ptr, 1, MPI_INT, 0,angle_comm );//hold the number of elements that each process will send the root process
            
            MPI_Barrier(angle_comm);
            MPI_Barrier(sub_world_comm);
            if (angle_id==0)
            {
                printf("Angle_procs: %d 1st gather: %d, %d, %d\n", angle_procs, *(total_cont_procs_angle_Ptr), *(total_cont_procs_angle_Ptr+1), *(total_cont_procs_angle_Ptr+2));
            }
            
            MPI_Reduce(&count_cont_procs, &total_cont_procs_angle, 1, MPI_INT, MPI_SUM, 0, angle_comm); //for each angle sum the number of procs to continue and pass it to the root for angle_comm
            
            cont_proc_ids_anglePtr=malloc(total_cont_procs_angle*sizeof(int)); //each root proc in angle comm has to hold the id's of the old set of processes to cont
            
            *(displPtr+0)=0;
            if (angle_id==0)
            {
                for (j=1;j<angle_procs;j++)
                {
                    *(displPtr+j)=(*(displPtr+j-1))+(*(total_cont_procs_angle_Ptr+j-1 )); //set the displacement for each proces to put its vector of pprocess IDs that need to be continued
                    printf("Displacement: %d\n", *(displPtr+j));
                }
                
            }
            
            MPI_Gatherv(cont_proc_idsPtr,count_cont_procs,MPI_INT, cont_proc_ids_anglePtr, total_cont_procs_angle_Ptr, displPtr , MPI_INT, 0,angle_comm ); //send the vectors with the ids of the old processes that need to be cont to root in angle_comm
            
            MPI_Barrier(angle_comm);
            MPI_Barrier(sub_world_comm);
            
            if (angle_id==0)
            {
                printf("Total Cont Procs: %d\n", total_cont_procs_angle);
                for (j=0;j<total_cont_procs_angle;j++)
                {
                    {
                        printf("Number: %d ID: %d\n",j,  *(cont_proc_ids_anglePtr+j));
                    }
                }
            }
            
            //each root for angle_comm has the number of processes each angle range needs to restart and the array of what the IDs of those processes used to be
            //now have to combine all that info for rank 0 in MPI_COMM_WORLD and then end it to all processes in MPI_COMM_WORLD
            //if (myid==0)
            {
                free(displPtr);
                displPtr=NULL;
                //initalize variables to hold all data
                
                
                each_num_to_restart_per_anglePtr=malloc(num_angles*sizeof(int));
                displPtr=malloc(num_angles*sizeof(int));
                *(displPtr+0)=0;
            }
            
            MPI_Barrier(angle_comm);
            MPI_Barrier(sub_world_comm);
            
            if (angle_id==0)
            {
                //this is the part where all the root processes of angle_comm transfer thier info to the root proc of MPI_WORLD
                MPI_Reduce(&total_cont_procs_angle, &total_num_to_restart, 1, MPI_INT, MPI_SUM, 0, root_angle_comm); //for each angle sum the number of procs to continue and pass it to the root for MPI_COMM_WORLD
            
                MPI_Gather(&total_cont_procs_angle,1,MPI_INT, each_num_to_restart_per_anglePtr, 1, MPI_INT, 0,root_angle_comm );//hold the number of elements that each process sent the root  for MPI_COMM_WORLD
            
            
                if (myid==0)
                {
                    for (j=1;j<num_angles;j++)
                    {
                        *(displPtr+j)=(*(displPtr+j-1))+(*(each_num_to_restart_per_anglePtr+j-1 )); //set the displacement for each proces to put its vector of pprocess IDs that need to be continued
                    }
                }
            
            
                all_cont_process_idPtr=malloc(total_num_to_restart*sizeof(int));
            
            
                MPI_Gatherv(cont_proc_ids_anglePtr, total_cont_procs_angle, MPI_INT, all_cont_process_idPtr, each_num_to_restart_per_anglePtr,   displPtr, MPI_INT, 0, root_angle_comm);
            }
            
            MPI_Barrier(angle_comm);
            MPI_Barrier(sub_world_comm);
            
            if (myid==0)
            {
                printf("Global Cont Procs: %d\n", total_num_to_restart);
                for (j=0;j<total_num_to_restart;j++)
                {
                    {
                        printf("Global ID: %d\n", *(all_cont_process_idPtr+j));
                    }
                }
            }
            
            //destroy the old comms
            MPI_Barrier(angle_comm);
            MPI_Barrier(sub_world_comm);
            //destroy current angle comm and recreate a new one 
            MPI_Comm_free(&root_angle_comm);
            MPI_Comm_free(&angle_comm);
            MPI_Comm_free(&sub_world_comm);
            MPI_Group_free(&sub_world_group);
            MPI_Group_free(&world_group);
            free(cont_proc_idsPtr);
            free(cont_proc_ids_anglePtr);
            free(total_cont_procs_angle_Ptr);
            free(displPtr);
            //free(each_num_to_restart_per_anglePtr);
            //free(all_cont_process_idPtr);
        }
        
        //send all of myid==0 data to all processes in MPI_COMM_WORLD 
        MPI_Bcast( &total_num_to_restart, 1, MPI_INT, 0, MPI_COMM_WORLD );
        
        if (total_num_to_restart>0)
        {
            if (myid != 0 )
            {
                printf("Proc: %d, Global Cont Procs: %d\n", myid, total_num_to_restart);
                //allocate data of appropriate size for all processes to hold the data from MPI_Bcast
                tmp=realloc(all_cont_process_idPtr,total_num_to_restart *sizeof(int));
                if (tmp!=NULL)
                {
                    all_cont_process_idPtr=tmp;
                }
                else
                {
                    printf("Error with reserving space to hold data about restarting process ID's\n");
                }
                //free(tmp);
                printf("Proc: %d, Num_angles: %d\n", myid, num_angles);
                tmp=realloc(each_num_to_restart_per_anglePtr, num_angles*sizeof(int));
                if (tmp!=NULL)
                {
                    each_num_to_restart_per_anglePtr=tmp;
                }
                else
                {
                    printf("Error with reserving space to hold data about restarting process numbers for each angle range\n");
                }
                //free(tmp);
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast( all_cont_process_idPtr, total_num_to_restart, MPI_INT, 0, MPI_COMM_WORLD );
            MPI_Bcast( each_num_to_restart_per_anglePtr, num_angles, MPI_INT, 0, MPI_COMM_WORLD );
            MPI_Bcast( &old_num_angle_procs, 1, MPI_INT, 0, MPI_COMM_WORLD );
        
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid==numprocs-1)
            {
                printf("Number of processes: %d\n", old_num_angle_procs);
                printf("restarting process numbers for each angle range: %d, %d, %d\n", *(each_num_to_restart_per_anglePtr), *(each_num_to_restart_per_anglePtr+1), *(each_num_to_restart_per_anglePtr+2));
            }
        
            //assign proper number of processes to each angle range to con't sims and then reset angle_id to original value from when simulation was first started
            color=0; //by default all processes have this value
        
            count=0;
            for (j=0;j<num_angles;j++)
            {
                if (myid>=count   &&   myid<count+(*(each_num_to_restart_per_anglePtr+j)) )
                {
                    color=j;
                }   
                count+=(*(each_num_to_restart_per_anglePtr+j));
                printf("Myid: %d, Color: %d, Count %d, Num To Start Per Angle: %d\n", myid, color, count, (*(each_num_to_restart_per_anglePtr+j)));
            }
            
            if (count!=numprocs)
            {
                //if the number of processes needed to continue the simulation is different from the number of processes in the mpiexec call exit
                printf('The simulation needs %d processes to properly continue. The number of processes initialized was %d.\nThe program is now exiting to prevent data corruption\n.', count, numprocs);
                exit(2);
            }
        
        
            MPI_Comm_split(MPI_COMM_WORLD, color , myid, &angle_comm);
            MPI_Comm_rank(angle_comm, &angle_id);
            MPI_Comm_size(angle_comm, &angle_procs);
        
            printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myid, numprocs, angle_id, angle_procs);
        
            angle_procs=old_num_angle_procs;
        
            //reset the angle for each process
            theta_jmin_thread= (*(thread_theta+  color)) *(M_PI/180);
            theta_jmax_thread= theta_jmin_thread+(delta_theta*(M_PI/180));
                
            //reset the angle_id for each process
            count=0;
            for (j=0;j<color;j++)
            {
                count+=(*(each_num_to_restart_per_anglePtr+j));
            }
        
            angle_id=(*(all_cont_process_idPtr+count+angle_id));
        
            snprintf(mc_dir,sizeof(flash_prefix),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this
        }
        else
        {
            //if there are no more processes to continue just break up processes normally  so they read in checkpoint files of completed processes and jump to merging files
            MPI_Comm_split(MPI_COMM_WORLD, myid/procs_per_angle , myid, &angle_comm);
            MPI_Comm_rank(angle_comm, &angle_id);
            MPI_Comm_size(angle_comm, &angle_procs);
        }
        free(all_cont_process_idPtr);
        free(each_num_to_restart_per_anglePtr);
    }
      
    MPI_Barrier(MPI_COMM_WORLD);
    
    if ((theta_jmin_thread >= 0) &&  (theta_jmax_thread <= (2*M_PI/180) )) //if within small angle (0-2 degrees) use _small inj_radius and frm2 have to think about this for larger domains
    {
        inj_radius=inj_radius_small;
        frm2=frm2_small;
        frm0=frm0_small;
        ph_weight_suggest=ph_weight_small;
    }
    else
    {
        inj_radius=inj_radius_large;
        frm2=frm2_large;
        frm0=frm0_large;
        ph_weight_suggest=ph_weight_large;
    }
    
    //make vector to hold the frames we are injecting in, vector should have (frm2-frm0)/angle_procs slots, if fps is const
        proc_frame_size=ceil((frm2-frm0)/ (float) angle_procs);
        frame_array=malloc(((frm2-frm0)+1)*sizeof(int));
        
        for (j=0;j<((frm2-frm0)+1); j++)
        {
            *(frame_array+j)=frm0+j ;
            //printf("proc: %d frame: %d\n", angle_id, *(frame_array+j));
        }
       
            
        
            //set this now incase there is no checkpoint file, then this wont be overwritten and the corretc values will be passed even if the user decides to restart
            framestart=(*(frame_array +(angle_id*proc_frame_size)));
            scatt_framestart=framestart;
                 
            if (angle_id != (angle_procs-1)) 
            {
                frm2=(*(frame_array +((angle_id*proc_frame_size) + proc_frame_size-1) )); //section off blocks of the frame_array to give to each angle_id
            }
            else
            {
                frm2=(*(frame_array + (frm2-frm0) )); //if angle_id is last give it the last set, even if its uneven
            }
                
            
            if (restrt=='c')
            {
                printf(">> mc.py:  Reading checkpoint\n");
                //#pragma omp critical
                
                    readCheckpoint(mc_dir, &phPtr, &frm2, &framestart, &scatt_framestart, &num_ph, &restrt, &time_now, angle_id, &angle_procs);
                
                /*
                for (i=0;i<num_ph;i++)
                {
                    printf("%e,%e,%e, %e,%e,%e, %e, %e\n",(phPtr+i)->p0, (phPtr+i)->p1, (phPtr+i)->p2, (phPtr+i)->p3, (phPtr+i)->r0, (phPtr+i)->r1, (phPtr+i)->r2, (phPtr+i)->num_scatt );
                }
                */
                if (restrt=='c')
                {
                    printf(">> Rank %d: Starting from photons injected at frame: %d out of %d\n", angle_id,framestart, frm2);
                    printf(">> Rank %d with angles %0.1lf-%0.1lf: Continuing scattering %d photons from frame: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,num_ph, scatt_framestart);
                    printf(">> Rank %d with angles %0.1lf-%0.1lf: The time now is: %e\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,time_now);
                }
                else
                {
                    printf(">> Rank %d with angles %0.1lf-%0.1lf: Continuing simulation by injecting photons at frame: %d out of %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,framestart, frm2); //starting with new photon injection is same as restarting sim
                }
                
            }
            else if ((stat(mc_dir, &st) == -1) && (restrt=='r'))
            {
                mkdir(mc_dir, 0777); //make the directory with full permissions
                
                
            }
            else 
            {
                if (angle_id==0)
                {
                    printf(">> proc %d with angles %0.1lf-%0.1lf:  Cleaning directory \n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
                    dirp = opendir(mc_dir);
                    while ((entry = readdir(dirp)) != NULL) 
                    {
                        if (entry->d_type == DT_REG) { /* If the entry is a regular file */
                             file_count++; //count how many files are in dorectory
                        }
                    }
                    printf("File count %d\n", file_count);
                    
                    if (file_count>0)
                    {
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mc_proc_*"); //prepares string to remove *.dat in mc_dir
                        system(mc_operation);
                        
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mcdata_PW_*"); //prepares string to remove *.dat in mc_dir
                        system(mc_operation);
                        
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mcdata_PW*"); //prepares string to remove *.dat in mc_dir
                        system(mc_operation);
                        
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mc_chkpt_*.dat"); //prepares string to remove *.dat in mc_dir
                        system(mc_operation);
                        
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mc_output_*.log"); //prepares string to remove *.log in mc_dir
                        system(mc_operation);
                        
                    }
                }
            }
    
            #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
            if (framestart>=3000)
            {
                increment_inj=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                fps_modified=1; //therefore dt between files become 1 second
                
            }
            #else
            {
                increment_inj=1;
                fps_modified=fps;
            }
            #endif
                        
            dt_max=1.0/fps_modified;
           
            MPI_Barrier(angle_comm); 
            snprintf(log_file,sizeof(log_file),"%s%s%d%s",mc_dir,"mc_output_", angle_id,".log" );
            printf("%s\n",log_file);
            fPtr=fopen(log_file, "a");
            
            printf( "Im Proc %d with angles %0.1lf-%0.1lf proc_frame_size is %d Starting on Frame: %d Injecting until %d scatt_framestart: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, proc_frame_size, framestart, frm2, scatt_framestart);
            
            fprintf(fPtr, "Im Proc %d with angles %0.1lf-%0.1lf  Starting on Frame: %d scatt_framestart: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, framestart, scatt_framestart);
            fflush(fPtr);
            
            free(frame_array);
            
            //for a checkpoint implementation, start from the last saved "frame" value and go to the saved "frm2" value
            
            //#pragma omp for 
            
            for (frame=framestart;frame<=frm2;frame=frame+increment_inj)
            {
                
                #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
                if (frame>=3000)
                {
                    increment_inj=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                    fps_modified=1; //therefore dt between files become 1 second
                    
                }
                #else
                {
                    increment_inj=1;
                    fps_modified=fps;
                }
                #endif
                                
                 if (restrt=='r')
                 {
                    time_now=frame/fps; //for a checkpoint implmentation, load the saved "time_now" value when reading the ckeckpoint file otherwise calculate it normally
                 }
                
                //printf(">> mc.py: Working on Frame %d\n", frame);
                fprintf(fPtr,"Im Proc: %d with angles %0.1lf - %0.1lf Working on Frame: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, frame);
                fflush(fPtr);
                
                if (restrt=='r')
                {
                    
                    
                    //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                    #if DIMENSIONS == 2
                    {
                        //if (strcmp(flash_sim, this_sim)==0)
                        #if SIM_SWITCH == FLASH
                        //{
                            //if using FLASH data for 2D
                            //put proper number at the end of the flash file
                            modifyFlashName(flash_file, flash_prefix, frame);
                            
                            fprintf(fPtr,">> Im Proc: %d with angles %0.1lf-%0.1lf: Opening FLASH file %s\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, flash_file);
                            fflush(fPtr);
                            
                            readAndDecimate(flash_file, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                                    &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, min_theta, max_theta, fPtr);
                        //}
                        //else if (strcmp(pluto_amr_sim, this_sim)==0)
                        #elif SIM_SWITCH == PLUTO_CHOMBO
                        //{
                            modifyPlutoName(flash_file, flash_prefix, frame);
                            
                            fprintf(fPtr,">> Im Proc: %d with angles %0.1lf-%0.1lf: Opening PLUTO file %s\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, flash_file);
                            fflush(fPtr);
                            
                            readPlutoChombo(flash_file, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                                    &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, min_theta, max_theta, fPtr);
                            
                            //exit(0);
                        //}
                        #else
                        //{
                            //if using RIKEN hydro data for 2D szx becomes delta r szy becomes delta theta
                            readHydro2D(FILEPATH, frame, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                                        &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, fPtr);
                            //fprintf(fPtr, "%d\n\n", array_num);
                        //}
                        #endif
                        fprintf(fPtr, "Number of Hydro Elements %d\n", array_num);
        //exit(0);
                    }
                    #else
                    {
                        fprintf(fPtr,">> Im Proc: %d with angles %0.1lf-%0.1lf\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
                        fflush(fPtr);
                        
                        read_hydro(FILEPATH, frame, inj_radius, &xPtr,  &yPtr, &zPtr,  &szxPtr, &szyPtr, &rPtr,\
                                   &thetaPtr, &phiPtr, &velxPtr,  &velyPtr, &velzPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, fps_modified, fPtr);
                    }
                    #endif
                    
                    //check for run type
                    //if(strcmp(cyl, this_run)==0)
                    #if SIMULATION_TYPE == CYLINDRICAL_OUTFLOW
                    {
                        //printf("In cylindrical prep\n");
                        cylindricalPrep(gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num);
                    }
                    //else if (strcmp(sph, this_run)==0)
                    #elif SIMULATION_TYPE == SPHERICAL_OUTFLOW
                    {
                        //printf("In Spherical\n");
                        sphericalPrep(rPtr, xPtr, yPtr,gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num , fPtr);
                    }
                    //else if (strcmp(struct_sph, this_run)==0)
                    #elif SIMULATION_TYPE == STRUCTURED_SPHERICAL_OUTFLOW
                    {
                        //printf("In Structural Spherical\n");
                        structuredFireballPrep(rPtr, thetaPtr, xPtr, yPtr,gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num , fPtr);
                    }
                    #endif
                        
                    //determine where to place photons and how many should go in a given place
                    //for a checkpoint implmentation, dont need to inject photons, need to load photons' last saved data 
                    fprintf(fPtr,">>  Proc: %d with angles %0.1lf-%0.1lf: Injecting photons\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
                    fflush(fPtr);
                    
                    //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                    #if DIMENSIONS == 2
                    {
                        photonInjection(&phPtr, &num_ph, inj_radius, ph_weight_suggest, min_photons, max_photons,spect, array_num, fps_modified, theta_jmin_thread, theta_jmax_thread, xPtr, yPtr, szxPtr, szyPtr,rPtr,thetaPtr, tempPtr, velxPtr, velyPtr,rng, fPtr );
                    }
                    #else
                    {
                        photonInjection3D(&phPtr, &num_ph, inj_radius, ph_weight_suggest, min_photons, max_photons,spect, array_num, fps_modified, theta_jmin_thread, theta_jmax_thread, xPtr, yPtr, zPtr, szxPtr, szyPtr,rPtr,thetaPtr, phiPtr, tempPtr, velxPtr, velyPtr, velzPtr, rng, fPtr);

                    }
                    #endif
                    
                    //printf("This many Photons: %d\n",num_ph); //num_ph is one more photon than i actually have
                    
                    //for (i=0;i<num_ph;i++)
                    //    printf("%e,%e,%e \n",(phPtr+i)->r0, (phPtr+i)->r1, (phPtr+i)->r2 );
                    
                }
                
                all_time_steps=malloc(num_ph*sizeof(double));
                sorted_indexes=malloc(num_ph*sizeof(int));
                
                //scatter photons all the way thoughout the jet
                //for a checkpoint implmentation, start from the last saved "scatt_frame" value eh start_frame=frame or start_frame=cont_frame
                if (restrt=='r')
                {
                    scatt_framestart=frame; //have to make sure that once the inner loop is done and the outer loop is incrememnted by one the inner loop starts at that new value and not the one read by readCheckpoint()
                }
                
                for (scatt_frame=scatt_framestart;scatt_frame<=last_frm;scatt_frame=scatt_frame+increment_scatt)
                {
                    #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
                    if (scatt_frame>=3000)
                    {
                        increment_scatt=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                        fps_modified=1; //therefore dt between files become 1 second
                        
                    }
                    #else
                    {
                        increment_scatt=1;
                        fps_modified=fps;
                    }
                    #endif
                    
                    dt_max=1.0/fps_modified; //if working with RIKEN files and scatt_frame>=3000 dt  is 1 second between each subsequent frame
                    
                    fprintf(fPtr,">>\n");
                    fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Working on photons injected at frame: %d out of %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,frame, frm2);
                    fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Simulation type %d - Working on frame %d\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, SIMULATION_TYPE, scatt_frame);
                    fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Opening file...\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
                    fflush(fPtr);
                    
                    //set new seed to increase randomness?
                    gsl_rng_set(rng, gsl_rng_get(rng));
                    
                   
                    //if (strcmp(DIM_SWITCH, dim_2d_str)==0)
                    #if DIMENSIONS == 2
                    {
                        //if (strcmp(flash_sim, this_sim)==0)
                        #if SIM_SWITCH == FLASH
                        {
                            //if using FLASH data for 2D
                            //put proper number at the end of the flash file
                            modifyFlashName(flash_file, flash_prefix, scatt_frame);
                            
                            fprintf(fPtr,">> Im Proc: %d with angles %0.1lf-%0.1lf: Opening FLASH file %s\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, flash_file);
                            fflush(fPtr);
                            
                            readAndDecimate(flash_file, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                                    &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, min_theta, max_theta, fPtr);
                        }
                        //else if (strcmp(pluto_amr_sim, this_sim)==0)
                        #elif SIM_SWITCH == PLUTO_CHOMBO
                        {
                            modifyPlutoName(flash_file, flash_prefix, scatt_frame);
                            
                            fprintf(fPtr,">> Im Proc: %d with angles %0.1lf-%0.1lf: Opening PLUTO file %s\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, flash_file);
                            fflush(fPtr);
                            
                            readPlutoChombo(flash_file, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                                    &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, min_theta, max_theta, fPtr);
                            
                            exit(0);
                        }
                        #else
                        {
                            //if using RIKEN hydro data for 2D szx becomes delta r szy becomes delta theta
                            readHydro2D(FILEPATH, scatt_frame, inj_radius, fps_modified, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                                        &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, fPtr);
                            //fprintf(fPtr, "%d\n\n", array_num);
                        }
                        #endif
                        fprintf(fPtr, "Number of Hydo Elements %d\n", array_num);
        //exit(0);
                    }
                    #else
                    {
                        fprintf(fPtr,">> Im Proc: %d with angles %0.1lf-%0.1lf\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
                        fflush(fPtr);
                        
                        read_hydro(FILEPATH, frame, inj_radius, &xPtr,  &yPtr, &zPtr,  &szxPtr, &szyPtr, &rPtr,\
                                   &thetaPtr, &phiPtr, &velxPtr,  &velyPtr, &velzPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, 1, min_r, max_r, fps_modified, fPtr);
                    }
                    #endif
                    fprintf(fPtr, "Number of Hydo Elements %d\n", array_num);
                    
                    
                    //check for run type
                    //if(strcmp(cyl, this_run)==0)
                    #if SIMULATION_TYPE == CYLINDRICAL_OUTFLOW
                    {
                        //printf("In cylindrical prep\n");
                        cylindricalPrep(gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num);
                    }
                    //else if (strcmp(sph, this_run)==0)
                    #elif SIMULATION_TYPE == SPHERICAL_OUTFLOW
                    {
                        sphericalPrep(rPtr, xPtr, yPtr,gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num, fPtr );
                    }
                    //else if (strcmp(struct_sph, this_run)==0)
                    #elif SIMULATION_TYPE == STRUCTURED_SPHERICAL_OUTFLOW
                    {
                        //printf("In Structural Spherical\n");
                        structuredFireballPrep(rPtr, thetaPtr, xPtr, yPtr,gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num , fPtr);
                    }
                    #endif
                        //printf("The result of read and decimate are arrays with %d elements\n", array_num);
                        
                    fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: propagating and scattering %d photons\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,num_ph);
                    fflush(fPtr);
                    
                    frame_scatt_cnt=0;
                    find_nearest_grid_switch=1; // set to true so the function findNearestPropertiesAndMinMFP by default finds the index of the grid block closest to each photon since we just read in a file and the prior index is invalid
                    
                    while (time_now<((scatt_frame+increment_scatt)/fps))
                    {
                        //if simulation time is less than the simulation time of the next frame, keep scattering in this frame
                        //for RIKEN hydro data, theres still 10 fps but after frame 3000, file increment is 10 not 1, therefore modify dt_max not fps
                        
                        //go through each photon and find blocks closest to each photon and properties of those blocks to calulate mean free path
                        //and choose the photon with the smallest mfp and calculate the timestep
                        

                        ph_scatt_index=findNearestPropertiesAndMinMFP(phPtr, num_ph, array_num, hydro_domain_x, hydro_domain_y, &time_step, xPtr,  yPtr, zPtr, szxPtr, szyPtr, velxPtr,  velyPtr,  velzPtr, dens_labPtr, tempPtr,\
                                                                      all_time_steps, sorted_indexes, rng, find_nearest_grid_switch, fPtr);
                        
                        find_nearest_grid_switch=0; //set to zero (false) since we do not absolutely need to refind the index, this makes the function findNearestPropertiesAndMinMFP just check if the photon is w/in the given grid box still

                        
                        //fprintf(fPtr, "In main: %d, %e, Newest Method results: %d, %e\n", ph_scatt_index, time_step, *(sorted_indexes+0), *(all_time_steps+(*(sorted_indexes+0))) );
                        //fflush(fPtr);
                        //for (i=1;i<num_ph;i++)
                        //{
                        //    fprintf(fPtr, "Newest Method results: %d, %e\n", *(sorted_indexes+i), *(all_time_steps+(*(sorted_indexes+i))) );
                        //}
                        
                        
                         if (time_step<dt_max)
                        {
                            //scatter the photon
                            //fprintf(fPtr, "Passed Parameters: %e, %e, %e\n", (ph_vxPtr), (ph_vyPtr), (ph_tempPtr));

                            time_step=photonScatter( phPtr, num_ph, dt_max, all_time_steps, sorted_indexes, velxPtr, velyPtr,  velzPtr, tempPtr,  &ph_scatt_index, &frame_scatt_cnt, rng, fPtr );
                            time_now+=time_step;
                            
                            
                            
                            if (frame_scatt_cnt%1000 == 0)
                            {
                                fprintf(fPtr,"Scattering Number: %d\n", frame_scatt_cnt);
                                //fprintf(fPtr,"Scattering Photon Number: %d\n", ph_scatt_index);
                                fprintf(fPtr,"The local temp is: %e K\n", *(tempPtr + (phPtr+ph_scatt_index)->nearest_block_index) );
                                fprintf(fPtr,"Average photon energy is: %e ergs\n", averagePhotonEnergy(phPtr, num_ph)); //write function to average over the photons p0 can then do (1.6e-9) to get keV
                                fprintf(fPtr,"The last time step was: %e.\nThe time now is: %e\n", time_step,time_now);
                                fflush(fPtr);
                            }
                            //exit(0);
                        }
                        else
                        {
                            time_now+=dt_max;
                            
                            //for each photon update its position based on its momentum
                            
                            updatePhotonPosition(phPtr, num_ph, dt_max, fPtr);
                        }
                        
                        //printf("In main 2: %e, %d, %e, %e\n", ((phPtr+ph_scatt_index)->num_scatt), ph_scatt_index, time_step, time_now);

                    }
                    
                    //get scattering statistics
                    phScattStats(phPtr, num_ph, &max_scatt, &min_scatt, &avg_scatt, &avg_r);
                        
                    fprintf(fPtr,"The number of scatterings in this frame is: %d\n", frame_scatt_cnt);
                    fprintf(fPtr,"The last time step was: %e.\nThe time now is: %e\n", time_step,time_now);
                    fprintf(fPtr,"The maximum number of scatterings for a photon is: %d\nThe minimum number of scattering for a photon is: %d\n", max_scatt, min_scatt);
                    fprintf(fPtr,"The average number of scatterings thus far is: %lf\nThe average position of photons is %e\n", avg_scatt, avg_r);
                    fflush(fPtr);
                    
                    fprintf(fPtr, ">> Proc %d with angles %0.1lf-%0.1lf: Making checkpoint file\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
                    fflush(fPtr);
                
                    fprintf(fPtr, " mc_dir: %s\nframe %d\nfrm2: %d\nscatt_frame: %d\n num_photon: %d\ntime_now: %e\nlast_frame: %d\n", mc_dir, frame, frm2, scatt_frame, num_ph, time_now, last_frm  );
                    fflush(fPtr);

                    save_chkpt_success=saveCheckpoint(mc_dir, frame, frm2, scatt_frame, num_ph, time_now, phPtr, last_frm, angle_id, old_num_angle_procs);
                    
                    if (save_chkpt_success==0)
                    {
                        //if we saved the checkpoint successfully also save the photons to the hdf5 file, else there may be something wrong with the file system
                        printPhotons(phPtr, num_ph,  scatt_frame , frame, mc_dir, angle_id, fPtr);
                    }
                    else
                    {
                        fprintf(fPtr, "There is an issue with opening and saving the chkpt file therefore MCRaT is not saving data to the checkpoint or mc_proc files to prevent corruption of those data.\n");
                        printf("There is an issue with opening and saving the chkpt file therefore MCRaT is not saving data to the checkpoint or mc_proc files to prevent corruption of those data.\n");
                        fflush(fPtr);
                        exit(1);
                    }
                    //exit(0);
                    
                     //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
                    #if SIM_SWITCH == RIKEN && DIMENSIONS ==3
                    {
                        //if (RIKEN_SWITCH==1)
                        {
                            free(zPtr);free(phiPtr);free(velzPtr);
                            zPtr=NULL; phiPtr=NULL; velzPtr=NULL;
                        }
                    }
                    #endif
                
                    free(xPtr);free(yPtr);free(szxPtr);free(szyPtr);free(rPtr);free(thetaPtr);free(velxPtr);free(velyPtr);free(densPtr);free(presPtr);
                    free(gammaPtr);free(dens_labPtr);free(tempPtr);
                    xPtr=NULL; yPtr=NULL;  rPtr=NULL;thetaPtr=NULL;velxPtr=NULL;velyPtr=NULL;densPtr=NULL;presPtr=NULL;gammaPtr=NULL;dens_labPtr=NULL;
                    szxPtr=NULL; szyPtr=NULL; tempPtr=NULL;
                }
                
                restrt='r';//set this to make sure that the next iteration of propogating photons doesnt use the values from the last reading of the checkpoint file
                free(phPtr); 
                phPtr=NULL;
                free(all_time_steps);
                all_time_steps=NULL;
                free(sorted_indexes);
                sorted_indexes=NULL;
            } 
            save_chkpt_success=saveCheckpoint(mc_dir, frame, frm2, scatt_frame, 0, time_now, phPtr, last_frm, angle_id, old_num_angle_procs); //this is for processes using the old code that didnt restart efficiently
            fprintf(fPtr, "Process %d has completed the MC calculation.\n", angle_id);
            fflush(fPtr);
            
            //exit(0);
                
        MPI_Barrier(angle_comm);
        
        //merge files from each worker thread within a directory
        {
            
             increment_scatt=1;
             file_count=0;
             
             //count number of files
             for (i=frm0;i<=last_frm;i=i+increment_scatt)
             {
                 
                //if ((RIKEN_SWITCH==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && (i>=3000))
                #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
                if (i>=3000)
                {
                    increment_scatt=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                }
                #endif
                file_count++;
             }
             
             //holds number of files for each process to merge
             MPI_Comm_size(angle_comm, &angle_procs); //to get the proper number of processes within the group
             MPI_Comm_rank(angle_comm, &angle_id); //reset the value of angle_id to what it should actualy be to properly distribute files to merge
             
             proc_frame_size=floor(file_count/ (float) angle_procs);
             frame_array=malloc(file_count*sizeof(int));
             proc_frame_array=malloc(angle_procs*sizeof(int)); //sets index of each proceesed acquired value
             element_num=malloc(angle_procs*sizeof(int));
             
             for (i=0;i<angle_procs;i++)
             {
                *(proc_frame_array+i)=i*proc_frame_size;
                *(element_num+i)=1;
             }
             
             //make vector with the files in order to pass them to each of the processes
             increment_scatt=1;
             file_count=0;
             for (i=frm0;i<=last_frm;i=i+increment_scatt)
             {
                //if ((RIKEN_SWITCH==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && (i>=3000))
                #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
                if (i>=3000)
                {
                    increment_scatt=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                }
                #endif

                *(frame_array+file_count)=i ;
                file_count++;
                //printf("file_count: %d frame: %d\n",  file_count-1, *(frame_array+file_count-1));
             }
             //pass  first frame number that each rpocess should start to merge, can calulate the file it should merge until
             MPI_Scatterv(frame_array, element_num, proc_frame_array, MPI_INT, &frm0, 1, MPI_INT, 0, angle_comm);
             
             //fprintf(fPtr, "Value: last_frm: ,%d\n", file_count);
             //fflush(fPtr);
             
             //make sure all files get merged by giving the rest to the last process
             if (angle_id==angle_procs-1)
             {
                proc_frame_size=file_count-proc_frame_size*(angle_procs-1); //for last process take over the remaining number of files
             }
             //calculate what the last file the preocess should merge up to
             i=0;
             last_frm=frm0;
             while(i<proc_frame_size)
             {
                //if ((RIKEN_SWITCH==1) && (strcmp(DIM_SWITCH, dim_3d_str)==0) && (last_frm>=3000))
                #if SIM_SWITCH == RIKEN && DIMENSIONS == 3
                if (last_frm>=3000)
                {
                    increment_scatt=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                }
                #else
                {
                    increment_scatt=1;
                }
                #endif
             
                last_frm+=increment_scatt;
                i++;
             }
             
             
            //if (angle_id==0)
            {
                //fprintf(fPtr, ">> Proc %d with angles %0.1lf-%0.1lf: Merging Files from %d to %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, frm0, last_frm);
                fprintf(fPtr, ">> Proc %d with angles %0.1lf-%0.1lf: Merging Files from %d to %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, frm0, last_frm);
                fflush(fPtr);
                
                dirFileMerge(mc_dir, frm0, last_frm, old_num_angle_procs, angle_id, fPtr);
            }
        }
        
        fprintf(fPtr, "Process %d has completed merging files.\n", angle_id);
        fflush(fPtr);
            
    fclose(fPtr);
    gsl_rng_free (rng);
   	 
    MPI_Finalize();
    //free(rng);
    //free(thread_theta);
    
	return 0;    
}
