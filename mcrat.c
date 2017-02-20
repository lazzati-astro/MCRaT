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
#include <omp.h>


#define THISRUN "Cylindrical"
#define FILEPATH "/Volumes/DATA6TB/Collapsars/2D/HUGE_BOXES/CONSTANT/16TI/"
#define FILEROOT "rhd_jet_big_13_hdf5_plt_cnt_"
#define MC_PATH "CMC_16TI_CYLINDRICAL/"
//#define MC_PATH "MC_16OI/Single_Photon_Cy_mc_total/"
#define MCPAR "mc.par"

int main(int argc, char **argv)
{
    //compile each time a macro is changed, have to supply the subfolder within the MC_PATH directory as a command line argument to the C program eg. MCRAT 1/
    
	// Define variables
	char flash_prefix[200]="";
	char mc_file[200]="" ;
    char this_run[200]=THISRUN;
    char *cyl="Cylindrical";
    char *sph="Spherical";
    char spect;//type of spectrum
    char restrt;//restart or not
    double fps, theta_jmin, theta_jmax ;//frames per second of sim, min opening angle of jet, max opening angle of jet in radians
    double inj_radius_small, inj_radius_large,  ph_weight_suggest ;//radius at chich photons are injected into sim
    int frm0,last_frm, frm2_small, frm2_large, j=0, min_photons, max_photons ;//frame starting from, last frame of sim, frame of last injection
    
    int num_thread=0, angle_count=0;
    int half_threads=floor((num_thread/2));
    double *thread_theta=NULL; //saves ranges of thetas for each thread to go through
    double delta_theta=0;
   
    //new muliple threads injecting and propagating photons
    const gsl_rng_type *rng_t;
    gsl_rng **rng;
    gsl_rng_env_setup();
    rng_t = gsl_rng_ranlxs0;

    
    //want to break up simulation by angle and injection frame & have each thread save data in its own folder 
    //have each thread check if its directory is made and if its restarting (delete evrything) or if its continuing with a previous simulation
    //the angle and the injection frames will be the names of mc_dir, therefore read mc.par first in MC_XXX directory
    
    //make strings of proper directories etc.
	snprintf(flash_prefix,sizeof(flash_prefix),"%s%s",FILEPATH,FILEROOT );
	//snprintf(mc_dir,sizeof(flash_prefix),"%s%s",FILEPATH,MC_PATH);
    //snprintf(mc_file,sizeof(flash_prefix),"%s%s",mc_dir, MCPAR);
    snprintf(mc_file,sizeof(flash_prefix),"%s%s%s",FILEPATH, MC_PATH,MCPAR);

    
    //printf(">> mc.py:  Reading mc.par\n");
    
    readMcPar(mc_file, &fps, &theta_jmin, &theta_jmax,&inj_radius_small,&inj_radius_large, &frm0,&last_frm ,&frm2_small, &frm2_large, &ph_weight_suggest, &min_photons, &max_photons, &spect, &restrt, &num_thread); //thetas that comes out is in radians
    
    //printf("small: %d large: %d weights: %e, %c\n", min_photons, max_photons, ph_weight_suggest, restrt);
    if (num_thread==0)
    {
        num_thread=omp_get_max_threads(); //if user specifies 0 in mc.par, default to max number of threads possible
    }
    
    half_threads=floor((num_thread/2));
    thread_theta=malloc(((half_threads)+1)*sizeof(double) );
    
    
    rng = (gsl_rng **) malloc((num_thread ) * sizeof(gsl_rng *)); 
     rng[0] = gsl_rng_alloc (rng_t); //initalize first random number generator to seed the others with random numbers
     
    for(j=1;j<num_thread;j++)
    {
        rng[j] = gsl_rng_alloc (rng_t);
        gsl_rng_set(rng[j],gsl_rng_get(rng[0])); 
    }
    
    //divide up angles and frame injections among threads DONT WANT NUMBER OF THREADS TO BE ODD
    //assign ranges to array that hold them
    delta_theta=(theta_jmax-theta_jmin)/(half_threads);
    *(thread_theta+0)=theta_jmin;
    //printf("%e\n", *(thread_theta+0));
    
    for (j=1;j<(half_threads+1); j++)
    {
        *(thread_theta+j)=*(thread_theta+(j-1))+delta_theta;
        //printf("%e\n", *(thread_theta+j)*180/M_PI);
    }

    //printf("half_threads: %d\n", half_threads);
    //start parallel section with half the threads possible and assign each thread to its value of theta by rewriting its private value of theta min and max that was read from mcrat
    omp_set_nested(1); //allow for nested parallelization
    //omp_set_dynamic(0);
    //omp_set_num_threads(num_thread);
    //#pragma omp parallel firstprivate(theta_jmin_thread, theta_jmax_thread, mc_dir, phPtr, framestart, scatt_framestart, num_ph, restrt, time_now, st, dirp, file_count, mc_filename, mc_operation, dt_max)
    //printf("%d\n", omp_get_num_threads() );
    
    #pragma omp parallel for num_threads(half_threads) private(angle_count)
    for (angle_count=(int) (theta_jmin*(180/M_PI)); angle_count< (int) (theta_jmax*(180/M_PI))  ;angle_count++ )
    {
        
        //printf("%d\n", omp_get_num_threads() );
        double inj_radius;
        int frm2;
        char mc_filename[200]="";
        char mc_filename_2[200]="";
        char mc_operation[200]="";
        char mc_dir[200]="" ;
        int file_count = 0;
        DIR * dirp;
        struct dirent * entry;
        struct stat st = {0};
        double theta_jmin_thread=0, theta_jmax_thread=0;
        
                
        theta_jmin_thread= (angle_count)*(M_PI/180);//(*(thread_theta+omp_get_thread_num() ));
        theta_jmax_thread= (angle_count+1)*(M_PI/180);//(*(thread_theta+omp_get_thread_num()+1 ));
        //printf("Thread %d: %0.1lf, %0.1lf \n %d %d\n", omp_get_thread_num(),  theta_jmin_thread*180/M_PI,  theta_jmax_thread*180/M_PI, frm2_small, frm2_large );
        
        snprintf(mc_dir,sizeof(flash_prefix),"%s%s%0.1lf-%0.1lf/",FILEPATH,MC_PATH, theta_jmin_thread*180/M_PI, theta_jmax_thread*180/M_PI ); //have to add angle into this
        
        //printf(">> Thread %d in  MCRaT: I am working on path: %s \n",omp_get_thread_num(), mc_dir  );
        
        if ((theta_jmin_thread >= 0) &&  (theta_jmax_thread <= (2*M_PI/180) )) //if within small angle (0-2 degrees) use _small inj_radius and frm2
        {
            inj_radius=inj_radius_small;
            frm2=frm2_small;
        }
        else
        {
            inj_radius=inj_radius_large;
            frm2=frm2_large;
        }
        //printf("Thread %d: %0.1lf, %0.1lf \n %d %e\n", omp_get_thread_num(),  theta_jmin_thread*180/M_PI,  theta_jmax_thread*180/M_PI, frm2, inj_radius );

        //want to also have another set of threads that each has differing ranges of frame injections therefore do nested parallelism here so each thread can read its own checkpoint file
        //#pragma omp parallel num_threads(2) firstprivate(restrt)
        {
            char flash_file[200]="";
            char log_file[200]="";
            FILE *fPtr=NULL; //pointer to log file for each thread
            double *xPtr=NULL,  *yPtr=NULL,  *rPtr=NULL,  *thetaPtr=NULL,  *velxPtr=NULL,  *velyPtr=NULL,  *densPtr=NULL,  *presPtr=NULL,  *gammaPtr=NULL,  *dens_labPtr=NULL;
            double *szxPtr=NULL,*szyPtr=NULL, *tempPtr=NULL; //pointers to hold data from FLASH files
            int num_ph=0, array_num=0, ph_scatt_index=0, max_scatt=0, min_scatt=0,i=0; //number of photons produced in injection algorithm, number of array elleemnts from reading FLASH file, index of photon whch does scattering, generic counter
            double dt_max=0, thescatt=0, accum_time=0; 
            double  gamma_infinity=0, time_now=0, time_step=0, avg_scatt=0; //gamma_infinity not used?
            double ph_dens_labPtr=0, ph_vxPtr=0, ph_vyPtr=0, ph_tempPtr=0;// *ph_cosanglePtr=NULL ;
            int frame=0, scatt_frame=0, frame_scatt_cnt=0, scatt_framestart=0, framestart=0;
            struct photon *phPtr=NULL; //pointer to array of photons 
            
            
            //if (omp_get_thread_num()==0)
            //{
                //printf( "A:%d Im Thread: %d with ancestor %d working in %s\n", omp_get_num_threads(),  omp_get_thread_num(), omp_get_ancestor_thread_num(1), mc_dir);
            //}
            
            if (restrt=='c')
            {
                printf(">> mc.py:  Reading checkpoint\n");
                //#pragma omp critical
                {
                    readCheckpoint(mc_dir, &phPtr, frm0, &framestart, &scatt_framestart, &num_ph, &restrt, &time_now);
                
                /*
                for (i=0;i<num_ph;i++)
                {
                    printf("%e,%e,%e, %e,%e,%e, %e, %e\n",(phPtr+i)->p0, (phPtr+i)->p1, (phPtr+i)->p2, (phPtr+i)->p3, (phPtr+i)->r0, (phPtr+i)->r1, (phPtr+i)->r2, (phPtr+i)->num_scatt );
                }
                */
                if (restrt=='c')
                {
                    printf(">> Thread %d with ancestor %d: Starting from photons injected at frame: %d out of %d\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1),framestart, frm2);
                    printf(">> Thread %d with ancestor %d: Continuing scattering %d photons from frame: %d\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1),num_ph, scatt_framestart);
                    printf(">> Thread %d with ancestor %d: The time now is: %e\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1),time_now);
                }
                else
                {
                    printf(">> Thread %d with ancestor %d: Continuing simulation by injecting photons at frame: %d out of %d\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1),framestart, frm2); //starting with new photon injection is same as restarting sim
                }
                }
            }
            else if (stat(mc_dir, &st) == -1)
            {
                mkdir(mc_dir, 0777); //make the directory with full permissions
                
                framestart=frm0; //if restarting then start from parameters given in mc.par file
                scatt_framestart=frm0;
            }
            else 
            {
                //remove everything from MC directory to ensure no corruption of data if theres other files there besides the mc.par file
                //for a checkpoint implementation, need to find the latest checkpoint file and read it and not delete the files
                #pragma omp critical
                {
                    printf(">> Thread %d with ancestor %d:  Cleaning directory \n",omp_get_thread_num(), omp_get_ancestor_thread_num(1));
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
                        for (i=0;i<=last_frm;i++)
                        {
                            
                            snprintf(mc_filename,sizeof(mc_filename),"%s%s%d%s", mc_dir,"mcdata_",i,"_P0.dat");
                            if(( access( mc_filename, F_OK ) != -1 ) )
                            {
                                snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s%d%s","exec rm ", mc_dir,"mcdata_",i,"_*.dat"); //prepares string to remove *.dat in mc_dir
                                printf("%s\n",mc_operation);
                                system(mc_operation);
                                
                                //snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s%d%s","exec rm ", mc_dir,"mcdata_",i,"_*");
                                //system(mc_operation);
                            }
                        }
                        //snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mcdata_PW_*.dat"); //prepares string to remove *.dat in mc_dir
                        //system(mc_operation);
                        
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mcdata_PW.dat"); //prepares string to remove *.dat in mc_dir
                        system(mc_operation);
                        
                        //snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mc_chkpt_*.dat"); //prepares string to remove *.dat in mc_dir
                        //system(mc_operation);
                        
                        snprintf(mc_operation,sizeof(flash_prefix),"%s%s%s","exec rm ", mc_dir,"mc_output_*.log"); //prepares string to remove *.log in mc_dir
                        system(mc_operation);
                        
                    }
                }
                framestart=frm0; //if restarting then start from parameters given in mc.par file
                scatt_framestart=frm0;
                
            }
            
            dt_max=1.0/fps;
            //#pragma omp barrier
            
            snprintf(log_file,sizeof(log_file),"%s%s",mc_dir,"mc_output.log" );
            printf("%s\n",log_file);
            fPtr=fopen(log_file, "w");
            
            fprintf(fPtr, "%d Im Thread: %d with ancestor %d Starting on Frame: %d scatt_framestart: %d\n", omp_get_num_threads(),  omp_get_thread_num(), omp_get_ancestor_thread_num(1), framestart, scatt_framestart);
            fflush(fPtr);
            //fclose(fPtr);
            //#pragma omp barrier
            //exit(0);
            //loop over frames 
            //for a checkpoint implementation, start from the last saved "frame" value and go to the saved "frm2" value
            
            //#pragma omp for 
            for (frame=framestart;frame<=frm2;frame++)
            {
                 if (restrt=='r')
                 {
                    time_now=frame/fps; //for a checkpoint implmentation, load the saved "time_now" value when reading the ckeckpoint file otherwise calculate it normally
                 }
                
                //printf(">> mc.py: Working on Frame %d\n", frame);
                fprintf(fPtr,"%d Im Thread: %d with ancestor %d Working on Frame: %d\n", omp_get_num_threads(),  omp_get_thread_num(), omp_get_ancestor_thread_num(1), frame);
                fflush(fPtr);
                
                if (restrt=='r')
                {
                    //put proper number at the end of the flash file
                    if (frame<10)
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%.3d%d",flash_prefix,000,frame);
                    }
                    else if (frame<100)
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%.2d%d",flash_prefix,00,frame);
                    }
                    else if (frame<1000)
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%d%d",flash_prefix,0,frame);
                    }
                    else
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%d",flash_prefix,frame);
                    }
                    
                    fprintf(fPtr,">> Im Thread: %d with ancestor %d: Opening FLASH file %s\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1), flash_file);
                    fflush(fPtr);
                    
                    //exit(0);
                    //read in FLASH file
                    //for a checkpoint implmentation, dont need to read the file yet
                    #pragma omp critical
                    {
                        readAndDecimate(flash_file, inj_radius, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                            &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, fPtr);
                        
                    }
                    
                    //check for run type
                    if(strcmp(cyl, this_run)==0)
                    {
                        //printf("In cylindrical prep\n");
                        cylindricalPrep(gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num);
                    }
                    else if (strcmp(sph, this_run)==0)
                    {
                        sphericalPrep(rPtr, xPtr, yPtr,gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num );
                    }
                        
                    //determine where to place photons and how many should go in a given place
                    //for a checkpoint implmentation, dont need to inject photons, need to load photons' last saved data 
                    fprintf(fPtr,">>  Thread: %d with ancestor %d: Injecting photons\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1));
                    fflush(fPtr);
                    photonInjection(&phPtr, &num_ph, inj_radius, ph_weight_suggest, min_photons, max_photons,spect, array_num, fps, theta_jmin_thread, theta_jmax_thread, xPtr, yPtr, szxPtr, szyPtr,rPtr,thetaPtr, tempPtr, velxPtr, velyPtr,rng[omp_get_thread_num()] ); 
                    //printf("%d\n",num_ph); //num_ph is one more photon than i actually have
                    /*
                    for (i=0;i<num_ph;i++)
                        printf("%e,%e,%e \n",(phPtr+i)->r0, (phPtr+i)->r1, (phPtr+i)->r2 );
                    */
                }
                
                //scatter photons all the way thoughout the jet
                //for a checkpoint implmentation, start from the last saved "scatt_frame" value eh start_frame=frame or start_frame=cont_frame
                if (restrt=='r')
                {
                    scatt_framestart=frame; //have to make sure that once the inner loop is done and the outer loop is incrememnted by one the inner loop starts at that new value and not the one read by readCheckpoint()
                }
                
                for (scatt_frame=scatt_framestart;scatt_frame<=last_frm;scatt_frame++)
                {
                    fprintf(fPtr,">>\n");
                    fprintf(fPtr,">> Thread %d with ancestor %d : Working on photons injected at frame: %d out of %d\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1),frame, frm2);
                    fprintf(fPtr,">> Thread %d with ancestor %d: %s - Working on frame %d\n",omp_get_thread_num(), omp_get_ancestor_thread_num(1), THISRUN, scatt_frame);
                    fprintf(fPtr,">> Thread %d with ancestor %d: Opening file...\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1));
                    fflush(fPtr);
                    
                    
                    //put proper number at the end of the flash file
                    if (scatt_frame<10)
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%.3d%d",flash_prefix,000,scatt_frame);
                    }
                    else if (scatt_frame<100)
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%.2d%d",flash_prefix,00,scatt_frame);
                    }
                    else if (scatt_frame<1000)
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%d%d",flash_prefix,0,scatt_frame);
                    }
                    else
                    {
                        snprintf(flash_file,sizeof(flash_prefix), "%s%d",flash_prefix,scatt_frame);
                    }
                    
                    #pragma omp critical
                    {
                        readAndDecimate(flash_file, inj_radius, &xPtr,  &yPtr,  &szxPtr, &szyPtr, &rPtr,\
                            &thetaPtr, &velxPtr,  &velyPtr,  &densPtr,  &presPtr,  &gammaPtr,  &dens_labPtr, &tempPtr, &array_num, fPtr);
                    }
                    
                    //check for run type
                    if(strcmp(cyl, this_run)==0)
                    {
                        //printf("In cylindrical prep\n");
                        cylindricalPrep(gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num);
                    }
                    else if (strcmp(sph, this_run)==0)
                    {
                        sphericalPrep(rPtr, xPtr, yPtr,gammaPtr, velxPtr, velyPtr, densPtr, dens_labPtr, presPtr, tempPtr, array_num );
                    }
                        //printf("The result of read and decimate are arrays with %d elements\n", array_num);
                        
                    fprintf(fPtr,">> Thread %d with ancestor %d: propagating and scattering %d photons\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1),num_ph);
                    fflush(fPtr);
                    
                    frame_scatt_cnt=0;
                    while (time_now<((scatt_frame+1)/fps))
                    {
                        //if simulation time is less than the simulation time of the next frame, keep scattering in this frame
                        //go through each photon and find blocks closest to each photon and properties of those blocks to calulate mean free path
                        //and choose the photon with the smallest mfp and calculate the timestep
                        
                        ph_scatt_index=findNearestPropertiesAndMinMFP(phPtr, num_ph, array_num, &time_step, xPtr,  yPtr, velxPtr,  velyPtr,  dens_labPtr, tempPtr,\
                            &ph_dens_labPtr, &ph_vxPtr, &ph_vyPtr, &ph_tempPtr, rng[omp_get_thread_num()]);
                            
                        //printf("In main: %e, %d, %e, %e\n", *(ph_num_scatt+ph_scatt_index), ph_scatt_index, time_step, time_now);
                        //printf("In main: %e, %d, %e, %e\n",((phPtr+ph_scatt_index)->num_scatt), ph_scatt_index, time_step, time_now);
                        
                         if (time_step<dt_max)
                        {
                            
                            //update number of scatterings and time
                            //(*(ph_num_scatt+ph_scatt_index))+=1;
                            ((phPtr+ph_scatt_index)->num_scatt)+=1;
                            frame_scatt_cnt+=1;
                            time_now+=time_step;
                            
                            updatePhotonPosition(phPtr, num_ph, time_step);
                            
                            //scatter the photon
                            //printf("Passed Parameters: %e, %e, %e\n", (ph_vxPtr), (ph_vyPtr), (ph_tempPtr));

                            photonScatter( (phPtr+ph_scatt_index), (ph_vxPtr), (ph_vyPtr), (ph_tempPtr), rng[omp_get_thread_num()] );
                            
                            
                            if (frame_scatt_cnt%1000 == 0)
                            {
                                fprintf(fPtr,"Scattering Number: %d\n", frame_scatt_cnt);
                                fprintf(fPtr,"The local temp is: %e\n", (ph_tempPtr));
                                fprintf(fPtr,"Average photon energy is: %e\n", averagePhotonEnergy(phPtr, num_ph)); //write function to average over the photons p0 and then do (*3e10/1.6e-9)
                                fflush(fPtr);
                            }
                            
                        }
                        else
                        {
                            time_now+=dt_max;
                            
                            //for each photon update its position based on its momentum
                            
                            updatePhotonPosition(phPtr, num_ph, dt_max);
                        }
                        
                        //printf("In main 2: %e, %d, %e, %e\n", ((phPtr+ph_scatt_index)->num_scatt), ph_scatt_index, time_step, time_now);

                    }
                    
                //get scattering statistics
                phScattStats(phPtr, num_ph, &max_scatt, &min_scatt, &avg_scatt);
                        
                fprintf(fPtr,"The number of scatterings in this frame is: %d\n", frame_scatt_cnt);
                fprintf(fPtr,"The last time step was: %lf.\nThe time now is: %lf\n", time_step,time_now);
                fprintf(fPtr,"The maximum number of scatterings for a photon is: %d\nThe minimum number of scattering for a photon is: %d\n", max_scatt, min_scatt);
                fprintf(fPtr,"The average number of scatterings thus far is: %lf\n", avg_scatt);
                fflush(fPtr);
                
                printPhotons(phPtr, num_ph,  scatt_frame , frame, mc_dir);
                
                //for a checkpoint implmentation,save the checkpoint file here after every 5 frames or something
                //save the photons data, the scattering number data, the scatt_frame value, and the frame value
                //WHAT IF THE PROGRAM STOPS AFTER THE LAST SCATT_FRAME, DURING THE FIRST SCATT_FRAME OF NEW FRAME VARIABLE - save restrt variable as 'r'
                fprintf(fPtr, ">> Thread %d with ancestor %d: Making checkpoint file\n", omp_get_thread_num(), omp_get_ancestor_thread_num(1));
                fflush(fPtr);
                saveCheckpoint(mc_dir, frame, scatt_frame, num_ph, time_now, phPtr, last_frm);
                
                free(xPtr);free(yPtr);free(szxPtr);free(szyPtr);free(rPtr);free(thetaPtr);free(velxPtr);free(velyPtr);free(densPtr);free(presPtr);
                free(gammaPtr);free(dens_labPtr);free(tempPtr);
                xPtr=NULL; yPtr=NULL;  rPtr=NULL;thetaPtr=NULL;velxPtr=NULL;velyPtr=NULL;densPtr=NULL;presPtr=NULL;gammaPtr=NULL;dens_labPtr=NULL;
                szxPtr=NULL; szyPtr=NULL; tempPtr=NULL;
                }
                restrt='r';//set this to make sure that the next iteration of propogating photons doesnt use the values from the last reading of the checkpoint file
                free(phPtr); 
                phPtr=NULL;
                
            } 
        
        
        }//end omp parallel inner section
        
        //merge files from each worker thread within a directory
        //dirFileMerge(mc_dir, frm0, last_frm);
        
    } //end omp parallel section
    
    //gsl_rng_free (rand);
    //free rand number generator
    for (j=0;j<num_thread;j++)
    {
        gsl_rng_free(rng[j]);
    }
    free(rng);
    free(thread_theta);
    
	return 0;
    
}
