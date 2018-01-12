/*
  *  This code is to merge all the files among different directories once the MCRaT simulation is complete
  *  The input should be the main CMC directory and the sub directories of each angle range with the number of MPI processes used to start the simulation in each directory
  *  eg call: mpiexec -np X /.merge /dir/to/CMC_dir/  0.0-2.0/ 416 2.0-4.0/ 416 
  *  where X shuould be a multiple of the number of sub directories
  *  
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
//#include <gsl/gsl_rng.h>
//#include "mclib.h"
//#include "mpi.h"

void readMcPar(char file[200], double *fps, double *theta_jmin, double *theta_j, double *d_theta_j, double *inj_radius_small, double *inj_radius_large, int *frm0_small, int *frm0_large,\
int *last_frm, int *frm2_small,int *frm2_large, double *ph_weight_small,double *ph_weight_large,int *min_photons, int *max_photons, char *spect, char *restart, int *num_threads,  int *dim_switch);


#define MCPAR "mc.par"

int main(int argc, char **argv)
{
    int num_angle_dirs=0, i=0;
    int *num_procs_per_dir, frm0_small, frm0_large, last_frm, frm2_small, frm2_large, small_frm;
    int *frm_array;
    double garbage;
    char mc_file[200]="" ;
    
    char dir[300]="";
    char *str="mcdata_proc_";
    int  count=0;
    struct dirent* dent;
    int file_count = 0;
    DIR * dirp;
    struct dirent * entry;
    DIR* srcdir = opendir(argv[1]);

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
            num_angle_dirs++;
            printf("found directory %s\n", dent->d_name);
        }
        
    }
    
    closedir(srcdir);
    
    num_procs_per_dir=malloc(num_angle_dirs*sizeof(int));
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
            printf("found directory %s\n", dent->d_name);
            
            snprintf(dir,sizeof(dir),"%s%s/",argv[1],dent->d_name );
            dirs[count] =  malloc((strlen(dir)+1));
            strcpy(dirs[count],dir);
            printf("SECOND: found directory %s\n", dirs[count]);
            
            dirp = opendir(dir); //do into the directory to get each file
            while ((entry = readdir(dirp)) != NULL) 
            {
                if ((entry->d_type == DT_REG) && (strstr(entry->d_name, str) != NULL))
                { /* If the entry is a regular file  */
                    file_count++;
                    printf("%s\n", entry->d_name );
                }
            }
            *(num_procs_per_dir +count)=file_count;
            count++;
        }
        
    }
    
    closedir(srcdir);
    
    //find number of directories for each angle range
    printf("%s: %d\n", argv[1], num_angle_dirs);
    
    for (i=0;i<num_angle_dirs;i++)
    {
        printf("%d\n", *(num_procs_per_dir+i));
        printf(" %s\n",  dirs[i]);
    }
    
    //get the last and initial hydro file in sim
    snprintf(mc_file,sizeof(mc_file),"%s%s",argv[1],MCPAR);
    readMcPar(mc_file, &garbage,&garbage, &garbage, &garbage, &garbage,&garbage, &frm0_small,&frm0_large, &last_frm ,&frm2_small, &frm2_large, &garbage, &garbage, &i, &i, &i, &i, &i,&i); //thetas that comes out is in degrees
    
    printf("%s frm_0small: %d frm_0large: %d, last: %d\n", mc_file, frm0_small,frm0_large, last_frm);
    
    //with all the info make array of all the files that need to be created
    frm_array=malloc(sizeof(int)*(last_frm-small_frm));
    small_frm= (frm0_small < frm0_large) ? frm0_small : frm0_large;
    count=0;
    for (i=small_frm;i<last_frm+1;i++)
    {
        printf("Count: %d\n",i);
        *(frm_array+count)=i;
        count++;
    }
    
    //set up MPI and break up the processes into groups of the number of sub directories
    
    
}

void readMcPar(char file[200], double *fps, double *theta_jmin, double *theta_j, double *d_theta_j, double *inj_radius_small, double *inj_radius_large, int *frm0_small, int *frm0_large, int *last_frm, int *frm2_small,int *frm2_large , double *ph_weight_small,double *ph_weight_large,int *min_photons, int *max_photons, char *spect, char *restart, int *num_threads,  int *dim_switch)
{
    //function to read mc.par file
	FILE *fptr=NULL;
	char buf[100]="";
	double theta_deg;
	
	//open file
	fptr=fopen(file,"r");
	//read in frames per sec and other variables outlined in main()
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
    
    fscanf(fptr, "%d",num_threads);
    //printf("%d\n",*num_threads);
    fgets(buf, 100,fptr);
    
    fscanf(fptr, "%d",dim_switch);
    //printf("%d\n",*dim_switch);
    
	//close file
	fclose(fptr);
}



