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
 
* Version 9.0 late 2017 included full Klein Nishina Cross Section and polarization with stokes parameters
* Version 9.1 late 2018 including cyclosynchrotron absorption and emission
*/

#include "mcrat.h"

int main(int argc, char **argv)
{
    //compile each time a macro is changed, have to supply the subfolder within the MC_PATH directory as a command line argument to the C program eg. MCRAT 1/
    
	// Define variables
	char hydro_prefix[STR_BUFFER]="";
	char mc_file[STR_BUFFER]="" ;
    char spect;//type of spectrum
    char restrt;//restart or not
    double fps, fps_modified, theta_jmin, theta_jmax,hydro_domain_y, hydro_domain_x ;//frames per second of sim, min opening angle of jet, max opening angle of jet in radians, max y value in hydro domain
    double inj_radius_small, inj_radius_large,  ph_weight_suggest=1e50, ph_weight_small, ph_weight_large, *inj_radius_input=NULL, ph_weight_default=1e50;//radius at chich photons are injected into sim
    int frm0,last_frm, frm2_small, frm2_large, j=0, min_photons, max_photons, frm0_small, frm0_large, *frm2_input=NULL, *frm0_input=NULL ;//frame starting from, last frame of sim, frame of last injection
    int dim_switch=0;
    int find_nearest_grid_switch=0;
    int increment_inj=1, increment_scatt=1; //increments for injection loop and scattering loop, outer and inner loops respectively, the increment can change for RIKEN 3D hydro files
    
    double inj_radius;
    int frm2,save_chkpt_success=0;
    char mc_filename[STR_BUFFER]="";
    char mc_filename_2[STR_BUFFER]="";
    char mc_operation[STR_BUFFER]="";
    char mc_dir[STR_BUFFER]="" ;
    int file_count = 0;
    DIR * dirp;
    struct dirent * entry;
    struct stat st = {0};
    double theta_jmin_thread=0, theta_jmax_thread=0;
        
    char hydro_file[STR_BUFFER]="";
    char log_file[STR_BUFFER]="";
    FILE *fPtr=NULL; //pointer to log file for each thread
    double *xPtr=NULL,  *yPtr=NULL,  *rPtr=NULL,  *thetaPtr=NULL,  *velxPtr=NULL,  *velyPtr=NULL,  *densPtr=NULL,  *presPtr=NULL,  *gammaPtr=NULL,  *dens_labPtr=NULL;
    double *szxPtr=NULL,*szyPtr=NULL, *tempPtr=NULL; //pointers to hold data from FLASH files
    double *phiPtr=NULL, *velzPtr=NULL, *zPtr=NULL, *all_time_steps=NULL ;
    int num_ph=0, scatt_cyclosynch_num_ph=0, num_null_ph=0, array_num=0, ph_scatt_index=0, num_photons_find_new_element=0, max_scatt=0, min_scatt=0,i=0; //number of photons produced in injection algorithm, number of array elleemnts from reading FLASH file, index of photon whch does scattering, generic counter
    double dt_max=0, thescatt=0, accum_time=0; 
    double  gamma_infinity=0, time_now=0, time_step=0, avg_scatt=0,avg_r=0; //gamma_infinity not used?
    double ph_dens_labPtr=0, ph_vxPtr=0, ph_vyPtr=0, ph_tempPtr=0, ph_vzPtr=0;// *ph_cosanglePtr=NULL ;
    double min_r=0, max_r=0, min_theta=0, max_theta=0, nu_c_scatt=0, n_comptonized=0, remaining_time=0;
    int frame=0, scatt_frame=0, frame_scatt_cnt=0, frame_abs_cnt=0, scatt_framestart=0, framestart=0;
    struct photon *phPtr=NULL; //pointer to array of photons
    struct photon *scattered_photon=NULL; //pointer to the most recently scattered photon
    struct photonList photon_list; //pointer to array of photons
    struct hydro_dataframe hydrodata; //pointer to array of hydro data
    
    int angle_count=0, num_cyclosynch_ph_emit=0;
    int num_angles=0, old_num_angle_procs=0; //old_num_angle_procs is to hold the old number of procs in each angle when cont sims, if  restarting sims this gets set to angle_procs
    int *frame_array=NULL, *proc_frame_array=NULL, *element_num=NULL, *sorted_indexes=NULL,  proc_frame_size=0;
    double *thread_theta=NULL; //saves ranges of thetas for each thread to go through
    double delta_theta=1, num_theta_bins=0;
    double test_cyclosynch_inj_radius=0;
    
    int myid, numprocs, angle_procs, angle_id, procs_per_angle;
    int temporary[3]={0}, tempo=0;

    
   
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
    
    
    initalizePhotonList(&photon_list);
    hydroDataFrameInitialize(&hydrodata);
    
    readMcPar(&hydrodata, &theta_jmin, &theta_jmax, &num_theta_bins, &inj_radius_input, &frm0_input , &frm2_input, &min_photons, &max_photons, &spect, &restrt); //thetas that comes out is in degrees, need to free input frame and injection radius pointers
    fps=hydrodata.fps;//save this incase we need modifications to fps later on in hydro sim
    last_frm=hydrodata.last_frame;
    //printf("%c\n", restrt);
    
    //divide up angles and frame injections among threads DONT WANT NUMBER OF THREADS TO BE ODD
    //assign ranges to array that hold them
    
    //leave angles in degrees here
    num_angles= (int) num_theta_bins;
    delta_theta=((theta_jmax-theta_jmin)/num_theta_bins);
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
     if (restrt==INITALIZE) //uncomment this when I run MCRAT for sims that didnt originally save angle_procs
     {
        
        MPI_Comm_split(MPI_COMM_WORLD, myid/procs_per_angle , myid, &angle_comm);
        MPI_Comm_rank(angle_comm, &angle_id);
        MPI_Comm_size(angle_comm, &angle_procs);

        //printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myid, numprocs, angle_id, angle_procs);
                
        theta_jmin_thread= (*(thread_theta+  (myid/procs_per_angle))) *(M_PI/180);
        theta_jmax_thread= theta_jmin_thread+(delta_theta*(M_PI/180));

        snprintf(mc_dir,sizeof(mc_dir),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this

        old_num_angle_procs=angle_procs;
         
        inj_radius= (*(inj_radius_input+(myid/procs_per_angle)));
        frm2=(*(frm2_input+(myid/procs_per_angle)));
        frm0=(*(frm0_input+(myid/procs_per_angle)));
        ph_weight_suggest=ph_weight_default;

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
        
            snprintf(mc_dir,sizeof(mc_dir),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this
        
            //call the function to count the num of processes for each angle range that need to be con't
            int  count_cont_procs=0, total_cont_procs_angle=0, global_cont_procs=0;
            int *cont_proc_idsPtr=NULL, *total_cont_procs_angle_Ptr=NULL, *displPtr=NULL; //becomes the size of the number of old procceses 
            int *cont_proc_ids_anglePtr=NULL;
            
            old_num_angle_procs=getOrigNumProcesses(&count_cont_procs,  &cont_proc_idsPtr, mc_dir, angle_id,  angle_procs,  last_frm);
            
            //count_cont_procs=1;//just for testing purposes
            
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
            
            //if (angle_id==0)
            //{
            //    printf("Angle_procs: %d 1st gather: %d, %d, %d\n", angle_procs, *(total_cont_procs_angle_Ptr), *(total_cont_procs_angle_Ptr+1), *(total_cont_procs_angle_Ptr+2));
            //}
            
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
                printf("Proc %d: The total number of processes that still have work to do is: %d\n", myid, total_num_to_restart);
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
            //if (myid==numprocs-1)
            //{
            //    printf("Number of processes: %d\n", old_num_angle_procs);
                //printf("restarting process numbers for each angle range: %d, %d, %d\n", *(each_num_to_restart_per_anglePtr), *(each_num_to_restart_per_anglePtr+1), *(each_num_to_restart_per_anglePtr+2));
            //}
        
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
                printf("The simulation needs %d processes to properly continue. The number of processes initialized was %d.\nThe program is now exiting to prevent data corruption\n.", count, numprocs);
                exit(2);
            }
        
        
            MPI_Comm_split(MPI_COMM_WORLD, color , myid, &angle_comm);
            MPI_Comm_rank(angle_comm, &angle_id);
            MPI_Comm_size(angle_comm, &angle_procs);
        
            //printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", myid, numprocs, angle_id, angle_procs);
        
            angle_procs=old_num_angle_procs;
        
            //reset the angle for each process
            theta_jmin_thread= (*(thread_theta+  color)) *(M_PI/180);
            theta_jmax_thread= theta_jmin_thread+(delta_theta*(M_PI/180));
            
            inj_radius= (*(inj_radius_input+color));
            frm2=(*(frm2_input+color));
            frm0=(*(frm0_input+color));
            ph_weight_suggest=ph_weight_default;
                
            //reset the angle_id for each process
            count=0;
            for (j=0;j<color;j++)
            {
                count+=(*(each_num_to_restart_per_anglePtr+j));
            }
        
            angle_id=(*(all_cont_process_idPtr+count+angle_id));
        
            snprintf(mc_dir,sizeof(mc_dir),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this
        }
        else
        {
            //if there are no more processes to continue just break up processes normally  so they read in checkpoint files of completed processes and jump to merging files
            MPI_Comm_split(MPI_COMM_WORLD, myid/procs_per_angle , myid, &angle_comm);
            MPI_Comm_rank(angle_comm, &angle_id);
            MPI_Comm_size(angle_comm, &angle_procs);
            frm0=(*(frm0_input+0));
        }
        free(all_cont_process_idPtr);
        free(each_num_to_restart_per_anglePtr);
    }
      
    MPI_Barrier(MPI_COMM_WORLD);
    free(inj_radius_input); free(frm0_input);free(frm2_input); free(thread_theta); //free input arrays sonce we hav determined which of the values in the arrays are applicable for this process
    
    //make vector to hold the frames we are injecting in, vector should have (frm2-frm0)/angle_procs slots, if fps is const
    
    //angle_procs=1;//just for testing purposes

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
                
            
    if (restrt==CONTINUE)
    {
        printf(">> MCRaT: Reading checkpoint\n");
        //#pragma omp critical
        
            scatt_cyclosynch_num_ph=readCheckpoint(mc_dir, &photon_list, &frm2, &framestart, &scatt_framestart, &restrt, &time_now, angle_id, &angle_procs);
        
        /*
        for (i=0;i<num_ph;i++)
        {
            printf("%e,%e,%e, %e,%e,%e, %e, %e\n",(phPtr+i)->p0, (phPtr+i)->p1, (phPtr+i)->p2, (phPtr+i)->p3, (phPtr+i)->r0, (phPtr+i)->r1, (phPtr+i)->r2, (phPtr+i)->num_scatt );
        }
        */
        if (restrt==CONTINUE)
        {
            printf(">> Rank %d: Starting from photons injected at frame: %d out of %d\n", angle_id,framestart, frm2);
            printf(">> Rank %d with angles %0.1lf-%0.1lf: Continuing scattering %d photons from frame: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,photon_list.num_photons, scatt_framestart);
            printf(">> Rank %d with angles %0.1lf-%0.1lf: The time now is: %e\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,time_now);
        }
        else
        {
            printf(">> Rank %d with angles %0.1lf-%0.1lf: Continuing simulation by injecting photons at frame: %d out of %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,framestart, frm2); //starting with new photon injection is same as restarting sim
        }
        
    }
    else if ((stat(mc_dir, &st) == -1) && (restrt==INITALIZE))
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
            //file_count=0;
            
            if (file_count>0)
            {
                snprintf(mc_operation,sizeof(mc_operation),"%s%s%s","exec rm ", mc_dir,"mc_proc_*"); //prepares string to remove *.dat in mc_dir
                system(mc_operation);
                
                snprintf(mc_operation,sizeof(mc_operation),"%s%s%s","exec rm ", mc_dir,"mcdata_PW_*"); //prepares string to remove *.dat in mc_dir
                system(mc_operation);
                
                snprintf(mc_operation,sizeof(mc_operation),"%s%s%s","exec rm ", mc_dir,"mcdata_PW*"); //prepares string to remove *.dat in mc_dir
                system(mc_operation);
                
                snprintf(mc_operation,sizeof(mc_operation),"%s%s%s","exec rm ", mc_dir,"mc_chkpt_*.dat"); //prepares string to remove *.dat in mc_dir
                system(mc_operation);
                
                snprintf(mc_operation,sizeof(mc_operation),"%s%s%s","exec rm ", mc_dir,"mc_output_*.log"); //prepares string to remove *.log in mc_dir
                system(mc_operation);
                
            }
            
            free(dirp);
        }
    }
    
    #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
    if (framestart>=3000)
    {
        //increment_inj=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
        //fps_modified=1; //therefore dt between files become 1 second
        hydrodata.increment_inj_frame=10;
        hydrodata.fps=1;
    }
    #else
    {
        //increment_inj=1;
        //fps_modified=fps;
        hydrodata.increment_inj_frame=1;
        hydrodata.fps=fps; //this is already set in readMcPar function but may need to be modified
    }
    #endif
                
    dt_max=1.0/hydrodata.fps;
   
    MPI_Barrier(angle_comm);
    snprintf(log_file,sizeof(log_file),"%s%s%d%s",mc_dir,"mc_output_", angle_id,".log" );
    printf("%s\n",log_file);
    fPtr=fopen(log_file, "a");

    
    printf( "Im Proc %d with angles %0.1lf-%0.1lf proc_frame_size is %d Starting on Frame: %d Injecting until %d scatt_framestart: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, proc_frame_size, framestart, frm2, scatt_framestart);
    
    fprintf(fPtr, "Im Proc %d with angles %0.1lf-%0.1lf  Starting on Frame: %d scatt_framestart: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, framestart, scatt_framestart);
    fflush(fPtr);
    
    free(frame_array);
    frame_array=NULL;

    #if TAU_CALCULATION == TABLE
        fprintf(fPtr, "TAU_CALCULATION is set to TABLE\n");
        fflush(fPtr);

        //initalize the tabulated cross sections (if needed)
        initalizeHotCrossSection(myid, rng, fPtr);

        /*
        fprintf(fPtr, "Im Proc %d testing the hot cross section interpolation\n");
        fflush(fPtr);
        // Test interpolation (all ranks can do this now)
        double test1 = interpolateThermalHotCrossSection(log10(1e-2), 2.75, rng, fPtr);
        fprintf(fPtr, "Thermal test: %g %g %g\n", log10(1e-2), 2.75, test1);


        #if NONTHERMAL_E_DIST != OFF
            double test[N_GAMMA];
            interpolateSubgroupNonThermalHotCrossSection(log10(1e-2), test, rng, fPtr);
            fprintf(fPtr, "NonThermal test: %g %g %g %g\n", log10(1e-2), test[0], test[1], test[2]);
        #endif
        */
    #else
    fprintf(fPtr, "TAU_CALCULATION is set to DIRECT\n");
    fflush(fPtr);
    #endif

    
    //for a checkpoint implementation, start from the last saved "frame" value and go to the saved "frm2" value

    for (frame=framestart;frame<=frm2;frame=frame+hydrodata.increment_inj_frame)
    {
        hydrodata.inj_frame_number=frame;
        #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
        if (frame>=3000)
        {
            hydrodata.increment_inj_frame=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
            hydrodata.fps=1; //therefore dt between files become 1 second
            
        }
        #else
        {
            hydrodata.increment_inj_frame=1;
            hydrodata.fps=fps;
        }
        #endif
                        
         if (restrt==INITALIZE)
         {
            time_now=frame/hydrodata.fps; //for a checkpoint implmentation, load the saved "time_now" value when reading the ckeckpoint file otherwise calculate it normally
         }
        
        printHydroGeometry(fPtr);
        fprintf(fPtr,">> Im Proc: %d with angles %0.1lf - %0.1lf Working on Frame: %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, frame);
        fflush(fPtr);
        
        if (restrt==INITALIZE)
        {
            //can read FLASH 2D (no B field) and plutochombo and pluto dbl files in 2/2.5/3D with B field
            getHydroData(&hydrodata, frame, inj_radius, 1, min_r, max_r, min_theta, max_theta, fPtr);
                
            //determine where to place photons and how many should go in a given place
            //for a checkpoint implmentation, dont need to inject photons, need to load photons' last saved data
            fprintf(fPtr,">>  Proc: %d with angles %0.1lf-%0.1lf: Injecting photons\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
            fflush(fPtr);
            
            photonInjection(&photon_list, inj_radius, ph_weight_suggest, min_photons, max_photons,spect, theta_jmin_thread, theta_jmax_thread, &hydrodata,rng, fPtr );
            
            //printf("This many Photons: %d\n",num_ph); //num_ph is one more photon than i actually have
            
            //for (i=0;i<num_ph;i++)
            //    printf("%e,%e,%e \n",(phPtr+i)->r0, (phPtr+i)->r1, (phPtr+i)->r2 );
            
        }
        
        freeHydroDataFrame(&hydrodata);//free frame data here since we rewrite over pointers in next loop
        
        //scatter photons all the way thoughout the jet
        //for a checkpoint implmentation, start from the last saved "scatt_frame" value eh start_frame=frame or start_frame=cont_frame
        if (restrt==INITALIZE)
        {
            scatt_framestart=frame; //have to make sure that once the inner loop is done and the outer loop is incremented by one the inner loop starts at that new value and not the one read by readCheckpoint()
        }
        
        num_null_ph=0;
        hydrodata.increment_scatt_frame=1;
        for (scatt_frame=scatt_framestart;scatt_frame<=last_frm;scatt_frame=scatt_frame+hydrodata.increment_scatt_frame)
        {
            hydrodata.scatt_frame_number=scatt_frame;
            #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
            if (scatt_frame>=3000)
            {
                hydrodata.increment_scatt_frame=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
                hydrodata.fps=1; //therefore dt between files become 1 second
                
            }
            #else
            {
                hydrodata.increment_scatt_frame=1;
                hydrodata.fps=fps;
            }
            #endif
            
            dt_max=1.0/hydrodata.fps; //if working with RIKEN files and scatt_frame>=3000 dt  is 1 second between each subsequent frame
            
            fprintf(fPtr,">>\n");
            printHydroGeometry(fPtr);
            fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Working on photons injected at frame: %d out of %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,frame, frm2);
            
            #if SIMULATION_TYPE == SCIENCE
                fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Simulation type Science - Working on scattering photons in frame %d\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, scatt_frame);
            #elif SIMULATION_TYPE == SPHERICAL_OUTFLOW
                fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Simulation type Spherical Outflow - Working on scattering photons in frame %d\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, scatt_frame);
            #elif SIMULATION_TYPE == CYLINDRICAL_OUTFLOW
                fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Simulation type Cylindrical Outflow - Working on scattering photons in frame %d\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, scatt_frame);
            #elif SIMULATION_TYPE == STRUCTURED_SPHERICAL_OUTFLOW
                fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Simulation type Structured Spherical Outflow - Working on scattering photons in frame %d\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, scatt_frame);
            #endif
            
            //fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: Opening file...\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
            //fflush(fPtr);
            
            //set new seed to increase randomness?
            gsl_rng_set(rng, gsl_rng_get(rng));
            
            //calc min and max positions of photons
            phMinMax(&photon_list, &min_r, &max_r, &min_theta, &max_theta, fPtr);
            #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((scatt_frame != scatt_framestart) || (restrt==CONTINUE))
                //if ((scatt_frame == scatt_framestart) || (restrt==CONTINUE))//for testing
                {
                    //NEED TO DETERMINE IF min_r or max_r is smaller/larger than the rmin/rmax in photonEmitCyclosynch to properly emit photons in the range that the process is interested in
                    //printf("OLD: min_r %e max_r %e\n", min_r, max_r);
                    test_cyclosynch_inj_radius=calcCyclosynchRLimits( scatt_frame, frame, hydrodata.fps,  inj_radius, "min");
                    //printf("TEST MIN: %e\n", test);
                    min_r=(min_r < test_cyclosynch_inj_radius) ? min_r : test_cyclosynch_inj_radius ;
                    test_cyclosynch_inj_radius=calcCyclosynchRLimits( scatt_frame, frame, hydrodata.fps,  inj_radius, "max");
                    //printf("TEST MAX: %e\n", test);
                    max_r=(max_r > test_cyclosynch_inj_radius ) ? max_r : test_cyclosynch_inj_radius ;
                    //printf("NEW: min_r %e max_r %e\n", min_r, max_r);
                }
            #endif

            getHydroData(&hydrodata, scatt_frame, inj_radius, 0, min_r, max_r, min_theta, max_theta, fPtr);
            
            //emit synchrotron photons here
            num_cyclosynch_ph_emit=0;
            
            //by default want to allocat ememory for time_steps and sorted indexes to scatter
            // all_time_steps=malloc(num_ph*sizeof(double)); this is now in the photon struct
            //sorted_indexes=malloc(num_ph*sizeof(int)); this is no longer needed since creating the photonList struct
            
            #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((scatt_frame != scatt_framestart) || (restrt==CONTINUE)) //remember to revert back to !=
                //if ((scatt_frame == scatt_framestart) || (restrt==CONTINUE))//for testing
                {
                    //if injecting synch photons, emit them if continuing simulation from a point where scatt_frame != scatt_framestart
                    //if necessary, then add memory to then arrays allocated directly above
                    
                    fprintf(fPtr, "Emitting Cyclosynchrotron Photons in frame %d\n", scatt_frame);
                    
                    #if B_FIELD_CALC == INTERNAL_E
                        fprintf(fPtr, "Calculating the magnetic field using internal energy and epsilon_B is set to %lf.\n", EPSILON_B);
                    #elif B_FIELD_CALC == TOTAL_E
                        //otherwise calculate B from the total energy
                        fprintf(fPtr, "Calculating the magnetic field using the total energy and epsilon_B is set to %lf.\n", EPSILON_B);
                    #else
                        fprintf(fPtr, "Using the magnetic field from the hydro simulation.\n");
                    #endif
                    //fprintf(fPtr, "HYDRO_B_SCALE %lf.\n", HYDRO_B_SCALE);
                    
                    phScattStats(&photon_list, &max_scatt, &min_scatt, &avg_scatt, &avg_r, fPtr); //for testing synch photons being emitted where 'i' photons are

                    num_cyclosynch_ph_emit=photonEmitCyclosynch(&photon_list, inj_radius, ph_weight_suggest, max_photons, theta_jmin_thread, theta_jmax_thread, &hydrodata, rng, 0, 0, fPtr);
                }
            #endif
            
            fprintf(fPtr,">> Proc %d with angles %0.1lf-%0.1lf: propagating and scattering %d photons\n",angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI,photon_list.num_photons);
            fflush(fPtr);
            
            frame_scatt_cnt=0;
            frame_abs_cnt=0;
            find_nearest_grid_switch=1; // set to true so the function findNearestPropertiesAndMinMFP by default finds the index of the grid block closest to each photon since we just read in a file and the prior index is invalid
            num_photons_find_new_element=0;
            remaining_time=((scatt_frame+hydrodata.increment_scatt_frame)/hydrodata.fps)-time_now; //This keeps track of the amount of time that remains until a new frame has to be loaded, we initalize it to the time of the next frame and subtract the scattering times from it. If the time_now=time of the last hydro frame then this ==dt_max
            
            n_comptonized=0;
            //while (time_now<((scatt_frame+hydrodata.increment_scatt_frame)/hydrodata.fps))
            while (remaining_time>0)
            {
                //if simulation time is less than the simulation time of the next frame, keep scattering in this frame
                //for RIKEN hydro data, theres still 10 fps but after frame 3000, file increment is 10 not 1, therefore modify dt_max not fps
                
                //go through each photon and find blocks closest to each photon and properties of those blocks to calulate mean free path
                //and choose the photon with the smallest mfp and calculate the timestep
                num_photons_find_new_element+=findContainingHydroCell(&photon_list, &hydrodata, find_nearest_grid_switch, rng, fPtr);

                //now calculate all the mean free paths (and get a sorted index array)
                calcMeanFreePath(&photon_list, &hydrodata, rng, fPtr);

                find_nearest_grid_switch=0; //set to zero (false) since we do not absolutely need to refind the index, this makes the function findNearestPropertiesAndMinMFP just check if the photon is w/in the given grid box still
                
                //if (*(all_time_steps+(*(sorted_indexes+0)))<remaining_time)
                
                if (getPhoton(&photon_list, photon_list.sorted_indexes[0])->time_to_scatter < remaining_time)
                {
                    //scatter the photon
                    //fprintf(fPtr, "Passed Parameters: %e, %e, %e\n", (ph_vxPtr), (ph_vyPtr), (ph_tempPtr));

                    time_step=photonEvent( &photon_list, remaining_time, &hydrodata, &ph_scatt_index, &frame_scatt_cnt, &frame_abs_cnt, rng, fPtr );
                    time_now+=time_step;
                    
                    remaining_time-=time_step; //update the remaining time subtracting off the time that has been accumulated through photon scatterings.
                    
                    //get the photon so we can easily access of its properties
                    scattered_photon=getPhoton(&photon_list, ph_scatt_index);

                    
                    //see if the scattered phton was a seed photon, if so replenish the seed photon
                    #if CYCLOSYNCHROTRON_SWITCH == ON
                    if (scattered_photon->type == CS_POOL_PHOTON)
                    {
                        n_comptonized+=scattered_photon->weight;
                        scattered_photon->type = COMPTONIZED_PHOTON; //c for compton scattered synchrotron photon
                        
                        //fprintf(fPtr, "num_null_ph %d\n", num_null_ph);
                        //printf("The previous scattered photon was a seed photon %c.\n", scattered_photon->type);
                        num_cyclosynch_ph_emit+=photonEmitCyclosynch(&photon_list, inj_radius, ph_weight_suggest, max_photons, theta_jmin_thread, theta_jmax_thread, &hydrodata, rng, 1, ph_scatt_index, fPtr);
                        //fprintf(fPtr, " num_photon: %d\n",num_ph  );
                        //fflush(fPtr);
                        
                        scatt_cyclosynch_num_ph++;//keep track of the number of synch photons that have scattered for later in checking of we need to rebin them
                        //fprintf(fPtr,"photonEmitCyclosynch: scatt_cyclosynch_num_ph Number: %d\n", scatt_cyclosynch_num_ph);
                        //exit(0);
                         
                    }
                    #endif
                                        
                    if ((frame_scatt_cnt%1000 == 0) && (frame_scatt_cnt != 0)) //modified this so it doesn't print when all photons get absorbed at first and frame_scatt_cnt=0
                    {
                        fprintf(fPtr,"Scattering Number: %d\n", frame_scatt_cnt);
                        fprintf(fPtr,"The local temp is: %e K\n", *(hydrodata.temp + scattered_photon->nearest_block_index) );
                        fprintf(fPtr,"Average photon energy is: %e ergs\n", averagePhotonEnergy(&photon_list)); //write function to average over the photons p0 can then do (1.6e-9) to get keV
                        fprintf(fPtr,"The last time step was: %e.\nThe time now is: %e\n", time_step,time_now);
                        //fprintf(fPtr,"Before Rebin: The average number of scatterings thus far is: %lf\nThe average position of photons is %e\n", avg_scatt, avg_r);
                        fflush(fPtr);
                        
                        #if CYCLOSYNCHROTRON_SWITCH == ON
                        if (scatt_cyclosynch_num_ph>max_photons)
                        {
                            //if the number of synch photons that have been scattered is too high rebin them
                            
                            //printf("num_cyclosynch_ph_emit: %d\n", num_cyclosynch_ph_emit);
                            rebinCyclosynchCompPhotons(&photon_list, &num_cyclosynch_ph_emit, &scatt_cyclosynch_num_ph, max_photons, theta_jmin_thread, theta_jmax_thread, rng, fPtr);

                            //fprintf(fPtr, "rebinSynchCompPhotons: scatt_cyclosynch_num_ph: %d\n", scatt_cyclosynch_num_ph);
                            //exit(0);
                        }
                        #endif
                   }
                    //exit(0);
                }
                else
                {
                    //this handles the case of the scattering time being larger than the time to the next frame either right away or at any point in iteratively scattering the photons
                    time_now+=remaining_time;
                    
                    //for each photon update its position based on its momentum and the remaining time to the next frame
                    
                    updatePhotonPosition(&photon_list, remaining_time, fPtr);
                    
                    //if we are here, then we need to load the next frame and we can set the timestep=remaining_time and then set remaining_time=0
                    time_step=remaining_time;
                    remaining_time=0;
                }
                
                
                //printf("In main 2: %e, %d, %e, %e\n", (scattered_photon->num_scatt), ph_scatt_index, time_step, time_now);

            }
            
            #if CYCLOSYNCHROTRON_SWITCH == ON
            if ((scatt_frame != scatt_framestart) || (restrt==CONTINUE)) //rememebr to change to != also at the other place in the code
            //if ((scatt_frame == scatt_framestart) || (restrt==CONTINUE)) //for testing
            {
                if (scatt_cyclosynch_num_ph>max_photons)
                {
                    //rebin the photons to ensure that we have a constant amount here
                    fprintf(fPtr, "Num_ph: %d\n", photon_list.num_photons);
                    /*
                    fprintf(fPtr,"Before Rebin: The average number of scatterings thus far is: %lf\nThe average position of photons is %e\n", avg_scatt, avg_r);
                    fflush(fPtr);
                    */
                    rebinCyclosynchCompPhotons(&photon_list, &num_cyclosynch_ph_emit, &scatt_cyclosynch_num_ph, max_photons, theta_jmin_thread, theta_jmax_thread, rng, fPtr);
                  //exit(0);
               }
                                        

                
                //make sure the photons that shou;d be absorbed should be absorbed if we have actually emitted any synchrotron photons
                if (num_cyclosynch_ph_emit>0)
                {
                    n_comptonized-=phAbsCyclosynch(&photon_list, &frame_abs_cnt, &scatt_cyclosynch_num_ph, &hydrodata, fPtr);
                }
                
            }
            #endif
            
            //get scattering statistics
            phScattStats(&photon_list, &max_scatt, &min_scatt, &avg_scatt, &avg_r, fPtr);
                
            fprintf(fPtr,"The number of scatterings in this frame is: %d\n", frame_scatt_cnt);
            #if CYCLOSYNCHROTRON_SWITCH == ON
                fprintf(fPtr,"The number of cyclosynchrotron photons absorbed in this frame is: %d\n", frame_abs_cnt);
            #endif
            fprintf(fPtr,"The last time step was: %e.\nThe time now is: %e\n", time_step,time_now);
            fprintf(fPtr,"MCRaT had to refind the position of photons %d times in this frame.\n", num_photons_find_new_element);
            fprintf(fPtr,"The maximum number of scatterings for a photon is: %d\nThe minimum number of scatterings for a photon is: %d\n", max_scatt, min_scatt);
            fprintf(fPtr,"The average number of scatterings thus far is: %lf\nThe average position of photons is %e\n", avg_scatt, avg_r);
            
            fflush(fPtr);
            
            fprintf(fPtr, ">> Proc %d with angles %0.1lf-%0.1lf: Making checkpoint file\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI);
            fflush(fPtr);
        
            //fprintf(fPtr, " mc_dir: %s\nframe %d\nfrm2: %d\nscatt_frame: %d\n num_photon: %d\ntime_now: %e\nlast_frame: %d\n", mc_dir, frame, frm2, scatt_frame, num_ph, time_now, last_frm  );
            
            //fprintf(fPtr,"n_comptonized in this frame is: %e\n ", n_comptonized);
            //fflush(fPtr);
            
            save_chkpt_success=saveCheckpoint(mc_dir, frame, frm2, scatt_frame, time_now, &photon_list, last_frm, angle_id, old_num_angle_procs);
            
            if (save_chkpt_success==0)
            {
                //if we saved the checkpoint successfully also save the photons to the hdf5 file, else there may be something wrong with the file system
                printPhotons(&photon_list, frame_abs_cnt, num_cyclosynch_ph_emit, scatt_cyclosynch_num_ph, scatt_frame , frame, last_frm, mc_dir, angle_id, fPtr);
            }
            else
            {
                fprintf(fPtr, "There is an issue with opening and saving the chkpt file therefore MCRaT is not saving data to the checkpoint or mc_proc files to prevent corruption of those data.\n");
                printf("There is an issue with opening and saving the chkpt file therefore MCRaT is not saving data to the checkpoint or mc_proc files to prevent corruption of those data.\n");
                fflush(fPtr);
                exit(1);
            }
            
            #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
            {
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
            //free(all_time_steps); //malloc is called in for loop therefore free memory at th end of the loop
            //all_time_steps=NULL;
            free(sorted_indexes);
            sorted_indexes=NULL;
            freeHydroDataFrame(&hydrodata);
        }
        
        restrt=INITALIZE;//set this to make sure that the next iteration of propogating photons doesnt use the values from the last reading of the checkpoint file
        scatt_cyclosynch_num_ph=0; //set this back equal to 0 for next batch of injected/emitted photons starting from nect injection frame
        num_null_ph=0; //set this back equal to 0 for next batch of injected/emitted photons starting from nect injection frame
        free(phPtr);
        phPtr=NULL;
        freePhotonList(&photon_list);
        //free(all_time_steps);
        //all_time_steps=NULL;
        free(sorted_indexes);
        sorted_indexes=NULL;
    }
    save_chkpt_success=saveCheckpoint(mc_dir, frame, frm2, scatt_frame, time_now, &photon_list, last_frm, angle_id, old_num_angle_procs); //this is for processes using the old code that didnt restart efficiently

    fprintf(fPtr, "Process %d has completed the MC calculation.\n", angle_id);
    fflush(fPtr);
    
    //exit(0);
    #if TAU_CALCULATION == TABLE
        cleanupInterpolationData();
    #endif
                
    MPI_Barrier(angle_comm);
        
    //merge files from each worker thread within a directory

     hydrodata.increment_scatt_frame=1;
     file_count=0;
     
     //count number of files
     for (i=frm0;i<=last_frm;i=i+hydrodata.increment_scatt_frame)
     {
         
        #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
        if (i>=3000)
        {
            hydrodata.increment_scatt_frame=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
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
     hydrodata.increment_scatt_frame=1;
     file_count=0;
     for (i=frm0;i<=last_frm;i=i+hydrodata.increment_scatt_frame)
     {
        #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
        if (i>=3000)
        {
            hydrodata.increment_scatt_frame=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
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
        #if SIM_SWITCH == RIKEN && DIMENSIONS == THREE
        if (last_frm>=3000)
        {
            hydrodata.increment_scatt_frame=10; //when the frame ==3000 for RIKEN 3D hydro files, increment file numbers by 10 instead of by 1
        }
        #else
        {
            hydrodata.increment_scatt_frame=1;
        }
        #endif
     
        last_frm+=hydrodata.increment_scatt_frame;
        i++;
     }
     
         
    //fprintf(fPtr, ">> Proc %d with angles %0.1lf-%0.1lf: Merging Files from %d to %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, frm0, last_frm);
    fprintf(fPtr, ">> Proc %d with angles %0.1lf-%0.1lf: Merging Files from %d to %d\n", angle_id, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI, frm0, last_frm);
    fflush(fPtr);
    
    dirFileMerge(mc_dir, frm0, last_frm, old_num_angle_procs, angle_id, fPtr);

    fprintf(fPtr, "Process %d has completed merging files.\n", angle_id);
    fflush(fPtr);
            
    fclose(fPtr);
    gsl_rng_free (rng);
    
    free(frame_array);
    free(proc_frame_array);
    free(element_num);
    
    MPI_Finalize();
    
	return 0;    
}
