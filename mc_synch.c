/*
This file is for the different functions for emitting and absorbing synchrotron photons
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <glob.h>
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
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include "mclib.h"
#include <omp.h>
#include "mpi.h"

//#DEFINE CRITICAL_B FINE_STRUCT*sqrt(M_EL*C_LIGHT*C_LIGHT/pow(R_EL, 3))

double calcCyclotronFreq(double magnetic_field)
{
    //B has to be in gauss
    return CHARGE_EL*magnetic_field/(2*M_PI*M_EL*C_LIGHT);
}

double calcDimlessTheta(double temp)
{
    //temp has to be in kelvin
    return K_B*temp/(M_EL*C_LIGHT*C_LIGHT);
}

double calcB(double el_dens, double temp, double epsilon_b)
{
    //calc the B field from assuming its some fraction of the matter energy density
    return sqrt(8*M_PI*3*el_dens*K_B*temp/2);
}

double n_el_MJ(double el_dens, double dimlesstheta, double gamma)
{
    //function to calulate the number density of electrons using the maxwell juttner distribution 
    return el_dens*gamma*sqrt(gamma*gamma-1)*exp(-gamma/dimlesstheta)/(dimlesstheta*gsl_sf_bessel_Kn(2, 1.0/dimlesstheta));
}

double n_el_MB(double el_dens, double dimlesstheta, double gamma)
{
    //function to calc the number density of electrons at a given dimensionless temp and lorentz factor with the maxwell boltzmann dist
    double temp=dimlesstheta*(M_EL*C_LIGHT*C_LIGHT)/K_B;
    double v=C_LIGHT*sqrt(1-(1/pow(gamma, 2)));
    
    return el_dens*4*M_PI*pow(M_EL/(2*M_PI*K_B*temp) , 3/2)*(v*C_LIGHT*C_LIGHT/(pow(gamma, 3)))*exp((-M_EL*pow(v, 2))/(2*K_B*temp));
}

//These functions are to calculate the emissivity from Wardzinski+ 2000
double Z(double nu, double nu_c, double gamma )
{
    return pow(sqrt(pow(gamma,2)-1)*exp(1/gamma)/(1+gamma) ,2*nu*gamma/nu_c);
}

double Z_sec_der(double nu, double nu_c, double gamma)
{
    //calculated from mathematica and plugged in theta (from paper=pi/2)
    return nu*(-2*pow(gamma,3)*(1+gamma) + 4*pow(gamma,4)*(1+gamma-pow(gamma,2)-pow(gamma,3))*log(sqrt(pow(gamma,2)-1)*exp(1/gamma)/(1+gamma) ))/(nu_c*pow(gamma,5)*(1+gamma));
}

double chi(double dimlesstheta, double gamma)
{
    double val=0;
    
    if (dimlesstheta<=0.08)
    {
        val=sqrt(2*dimlesstheta*(pow(gamma,2)-1)/(gamma*(3*pow(gamma,2)-1)));
    }
    else
    {
        val=sqrt(2*dimlesstheta/(3*gamma));
    }
    
    return val;
}

double gamma0(double nu, double nu_c, double dimlesstheta)
{
    double val=0;

    if (dimlesstheta<=0.08)
    {
        val=sqrt(pow(1+(2*nu*dimlesstheta/nu_c)*(1+(9*nu*dimlesstheta/(2*nu_c))), (-1.0/3.0) ));
    }
    else
    {
        val=sqrt(pow((1+(4*nu*dimlesstheta/(3*nu_c))), (2.0/3.0)) );
    }
    
    return val;
}

double jnu(double nu, double nu_c, double dimlesstheta, double el_dens)
{
    double dimlesstheta_ref=calcDimlessTheta(1e7);
    double gamma=gamma0(nu, nu_c, dimlesstheta);
    double val=0;
    
    if (dimlesstheta<dimlesstheta_ref)
    {
        val=(pow(M_PI,(3.0/2.0))*pow(CHARGE_EL, 2)/(pow(2,(3.0/2.0))*C_LIGHT))*sqrt(nu*nu_c)*n_el_MB(el_dens, dimlesstheta, gamma)* Z(nu, nu_c, gamma)*chi( dimlesstheta, gamma)* pow(fabs(Z_sec_der(nu, nu_c, gamma)),(-1.0/2.0));
    }
    else
    {
        val=(pow(M_PI,(3.0/2.0))*pow(CHARGE_EL, 2)/(pow(2,(3.0/2.0))*C_LIGHT))*sqrt(nu*nu_c)*n_el_MJ(el_dens, dimlesstheta, gamma)* Z(nu, nu_c, gamma)*chi( dimlesstheta, gamma)* pow(fabs(Z_sec_der(nu, nu_c, gamma)),(-1.0/2.0));
    }
    
    return val;
}

double jnu_ph_spect(double nu, void *p)
{
    //struct jnu_params * params = (struct jnu_params *)p;
    double *params;
    params=(double *)p;
    double nu_c = params[0];//(params->nu_c);
    double dimlesstheta = params[1];//(params->dimlesstheta);
    double el_dens = params[2];//(params->el_dens);
    
    //printf("Nu %e nu_c %e dimlesstheta %e el_dens %e VAL %e\n", nu, nu_c, dimlesstheta, el_dens, jnu(nu, nu_c, dimlesstheta, el_dens)/(PL_CONST*nu));
    
    return jnu(nu, nu_c, dimlesstheta, el_dens)/(PL_CONST*nu);
}

double blackbody_ph_spect(double nu, void *p)
{
    //struct jnu_params * params = (struct jnu_params *)p;
    double *params;
    params=(double *)p;
    double temp = params[0];//(params->nu_c);
    
    //printf("Nu %e nu_c %e dimlesstheta %e el_dens %e VAL %e\n", nu, nu_c, dimlesstheta, el_dens, jnu(nu, nu_c, dimlesstheta, el_dens)/(PL_CONST*nu));
    
    return (PL_CONST*nu * pow((nu/C_LIGHT),2.0))/(exp(PL_CONST*nu/(K_B*temp))-1)/(PL_CONST*nu);
}

//the functions here are to calculate the total absorption cross section from Ghisellini+ 1991
double C(double nu_ph, double nu_c, double gamma_el, double p_el)
{
    return ((2.0*pow(gamma_el,2)-1)/(gamma_el*pow(p_el,2)))+2*nu_ph*((gamma_el/pow(p_el,2))-gamma_el*log((gamma_el+1)/p_el))/nu_c;
}

double G(double gamma_el, double p_el)
{
    return sqrt(1-2*pow(p_el,2)*(gamma_el*log((gamma_el+1)/p_el)-1));
}

double G_prime(double gamma_el, double p_el)
{
    return (3*gamma_el-(3*pow(gamma_el,2)-1)*log((gamma_el+1)/p_el))/G(gamma_el, p_el);
}

double synCrossSection(double el_dens, double T, double nu_ph, double p_el, double epsilon_b)
{
    double b_cr=FINE_STRUCT*sqrt(M_EL*C_LIGHT*C_LIGHT/pow(R_EL,3.0));
    double B=calcB(el_dens, T, epsilon_b);
    double nu_c=calcCyclotronFreq(B);
    double gamma_el=sqrt(p_el*p_el+1);
    
    //printf("calc gamma %e, temp %e\n", gamma_el, T);
    
    return (3.0*M_PI*M_PI/8.0)*(THOM_X_SECT/FINE_STRUCT)*(b_cr/B)*pow(nu_c/nu_ph, 2.0) * exp(-2*nu_ph*(gamma_el*log((gamma_el+1)/p_el)-1)/nu_c)* ((C(nu_ph, nu_c, gamma_el, p_el)/G(gamma_el, p_el))-(G_prime(gamma_el, p_el)/pow(G(gamma_el, p_el),2.0)));
}

double calcSynchRLimits(int frame_scatt, int frame_inj, double fps,  double r_inj, char *min_or_max)
{
    double val=r_inj;
    if (strcmp(min_or_max, "min")==0)
    {
        //printf("IN MIN\nframe_scatt %e frame_inj %e fps %e r_inj %e C_LIGHT %e\n", frame_scatt, frame_inj, fps, r_inj, C_LIGHT);
        val+= (C_LIGHT*(frame_scatt-frame_inj-1)/(2*fps));
    }
    else
    {
        //printf("IN MAX\n");
        val+= (C_LIGHT*(frame_scatt-frame_inj+1)/(2*fps));
    }
    
    //printf("Val %e\n", val);
    
    return val;
}

int photonEmitSynch(struct photon **ph_orig, int *num_ph, int *num_null_ph, double **all_time_steps, int **sorted_indexes, double r_inj, double ph_weight, int maximum_photons, int array_length, double fps, double theta_min, double theta_max , int frame_scatt, int frame_inj, double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temp, double *dens, double *vx, double *vy,  double epsilon_b, gsl_rng *rand, int riken_switch, int inject_single_switch, int scatt_ph_index, FILE *fPtr)
{
    double rmin=0, rmax=0, max_photons=0.1*maximum_photons; //have 10% as default, can change later need to figure out how many photons across simulations I want emitted
    double ph_weight_adjusted=0, position_phi=0;
    double dimlesstheta=0, nu_c=0, el_dens=0, error=0, ph_dens_calc=0, max_jnu=0;
    double el_p[4], ph_p_comv[4];
    double params[3];
    double fr_dum=0.0, y_dum=0.0, yfr_dum=0.0, com_v_phi=0, com_v_theta=0, position_rand=0;
    double *p_comv=NULL, *boost=NULL, *l_boost=NULL; //pointers to hold comov 4 monetum, the fluid vlocity, and the photon 4 momentum in the lab frame
    int status;
    int block_cnt=0, i=0, j=0, k=0, null_ph_count=0, *ph_dens=NULL, ph_tot=0, net_ph=0, min_photons=1;
    int *null_ph_indexes=NULL;
    int num_thread=omp_get_num_threads(), count_null_indexes=0, idx=0;
    struct photon *ph_emit=NULL; //pointer to array of structs that will hold emitted photon info
    struct photon *tmp=NULL;
    double *tmp_double=NULL;
    int *tmp_int=NULL;
    
    printf("IN EMIT SYNCH FUNCTION\n");
    fprintf(fPtr, "BEFORE Original number of photons: %d Null photons %d\n", (*num_ph), null_ph_count, ph_tot);
    fflush(fPtr);
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (10000);
    
    gsl_function F;
    F.function = &blackbody_ph_spect; //&jnu_ph_spect;
    
    
    if (inject_single_switch == 0)
    {
        rmin=calcSynchRLimits( frame_scatt, frame_inj, fps,  r_inj, "min");
        rmax=calcSynchRLimits( frame_scatt, frame_inj, fps,  r_inj, "max");
        
        //printf("rmin %e rmax %e, theta min/max: %e %e\n", rmin, rmax, theta_min, theta_max);
        #pragma omp parallel for num_threads(num_thread) reduction(+:block_cnt)
        for(i=0;i<array_length;i++)
        {
            
            //look at all boxes in width delta r=c/fps and within angles we are interested in NEED TO IMPLEMENT
            if ((*(r+i) >= rmin)  &&   (*(r+i)  < rmax  ) && (*(theta+i)< theta_max) && (*(theta+i) >=theta_min) )
            {
                block_cnt+=1;
            }
        }
        
        //fprintf(fPtr, "Block cnt %d\n", block_cnt);
        
        //allocate memory to record density of photons for each block
        ph_dens=malloc(block_cnt * sizeof(int));
        
        
        //calculate the photon density for each block and save it to the array
        j=0;
        ph_tot=-1;
        ph_weight_adjusted=ph_weight;
        while ((ph_tot>max_photons) || (ph_tot<min_photons) ) //can have 0 photons emitted
        {
            j=0;
            ph_tot=0;
            
            for (i=0;i<array_length;i++)
            {
                //printf("%d\n",i);
                //printf("%e, %e, %e, %e, %e, %e\n", *(r+i),(r_inj - C_LIGHT/fps), (r_inj + C_LIGHT/fps), *(theta+i) , theta_max, theta_min);
                if ((*(r+i) >= rmin)  &&   (*(r+i)  < rmax  ) && (*(theta+i)< theta_max) && (*(theta+i) >=theta_min) )
                {
                    if (riken_switch==0)
                    {
                        //using FLASH
                        //set parameters fro integration fo phtoons spectrum
                        el_dens= (*(dens+i))/M_P;
                        nu_c=calcCyclotronFreq(calcB(el_dens,*(temp+i) , epsilon_b));
                        dimlesstheta=calcDimlessTheta( *(temp+i));
                        //printf("Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e\n",*(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta);

                        params[0] = *(temp+i); //nu_c;
                        params[1]=dimlesstheta;
                        params[2]= el_dens;
                        F.params = &params;
                        
                        //printf("Integrating\n");
                        //status=gsl_integration_qags(&F, nu_c*1e-4, nu_c*1e2, 0, 1e-2, 10000, w, &ph_dens_calc, &error); //do this range b/c its ~4 order of magnitude difference between peak and lower limit and at high frequencies have exponential cut off therefore probability goes down very fast
                        //ph_dens_calc*=2*M_PI*(*(x+i))*pow(*(szx+i),2.0)/(fps*ph_weight_adjusted);
                        //printf ("error: %s\n", gsl_strerror (status));
                        //Before integrated the photon number spectrum to get the number of photons to emit
                        
                        //ph_p_comv[0]=nu_c*PL_CONST/C_LIGHT;
                        //ph_p_comv[1]=0;
                        //ph_p_comv[2]=0;
                        //ph_p_comv[3]=0; //dont care about these components of the 4 momentum
                        
                        //now assume steady state for thermal synchrotron radiation, emit number of photons as source function (j/absorption) at nu_c divided by h*nu_c
                        //singleElectron(&el_p[0], *(temp+i), &ph_p_comv[0], rand, fPtr); //get random electron, only care about its energy
                        //printf("Chosen el: p0 %e p1 %e p2 %e p3 %e\n", *(el_p+0), *(el_p+1), *(el_p+2), *(el_p+3));
                        //ph_dens_calc=jnu(nu_c, nu_c, dimlesstheta, el_dens)/(el_dens*synCrossSection(el_dens, *(temp+i), ph_p_comv[0]*C_LIGHT/PL_CONST, sqrt((el_p[0]*el_p[0]/(M_EL*M_EL*C_LIGHT*C_LIGHT))-1), epsilon_b)*(PL_CONST*nu_c));
                        
                        //printf("Integrating\n"); //instead integrating from 0 to nu_c
                        status=gsl_integration_qags(&F, 10, nu_c, 0, 1e-2, 10000, w, &ph_dens_calc, &error); //find number of low energy seed photons in the tail of the BB distribution
                        //printf ("error: %s\n", gsl_strerror (status));
                        
                        //printf("ph_dens_calc init=%e\n", ph_dens_calc);
                        ph_dens_calc*=2*M_PI*(*(x+i))*pow(*(szx+i),2.0)/(fps*ph_weight_adjusted);
                        
                        //printf("Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e, number of photons to emit %e, error %e, Intervals %zu\n",  *(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta, ph_dens_calc, error, w->size);
                        //exit(0);
                    }
                    else
                    {
                        printf("Emitting Photons with thermal synchrotron isn't available for non-FLASH non-2D hydro simulations.\n");
                        exit(0);
                    }
                    
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
            
            printf("dens: %d, photons: %d, adjusted weight: %e\n", *(ph_dens+(j-1)), ph_tot, ph_weight_adjusted);
            
        }
        
    }
    else
    {
        //do this when were emitting a synch photon to replace a scattered synchrotron photon
        ph_tot=1;
    }
    
    //FIND OUT WHICH PHOTONS IN ARRAY ARE OLD/WERE ABSORBED AND IDENTIFY THIER INDEXES AND HOW MANY, dont subtract this from ph_tot @ the end, WILL NEED FOR PRINT PHOTONS
    #pragma omp parallel for num_threads(num_thread) reduction(+:null_ph_count)
    for (i=0;i<*num_ph;i++)
    {
        if ((*ph_orig)[i].weight == 0)
        {
            null_ph_count+=1;
        }
    }
    
    
    fprintf(fPtr, "Original number of photons: %d Null photons %d\nEmitting %d synchrotron photons between %e and %e in this frame\n", (*num_ph), null_ph_count, ph_tot, rmin, rmax);
    fflush(fPtr);
    //exit(0);
    //allocate memory for that many photons and also allocate memory to hold comoving 4 momentum of each photon and the velocity of the fluid
    //ph_emit=malloc (ph_tot * sizeof (struct photon ));
    p_comv=malloc(4*sizeof(double));
    boost=malloc(4*sizeof(double));
    l_boost=malloc(4*sizeof(double));
    
    if (null_ph_count < ph_tot)
    {
        //if the totoal number of photons to be emitted is larger than the number of null phtons curently in the array, then have to grow the array
        //need to realloc memory to hold the old photon info and the new emitted photon's info
        fprintf(fPtr, "Allocating %d space\n", ((*num_ph)+ph_tot-null_ph_count));
        tmp=realloc(*ph_orig, ((*num_ph)+ph_tot-null_ph_count)* sizeof (struct photon )); //may have to look into directly doubling (or *1.5) number of photons each time we need to allocate more memory, can do after looking at profiling for "just enough" memory method
        if (tmp != NULL)
        {
            /* everything ok */
            *ph_orig = tmp;
            /*
             for (i=0;i<*num_ph;i++)
             {
             fprintf(fPtr, "i: %d after realloc freq %e\n", i, (*ph_orig)[i].p0*C_LIGHT/PL_CONST );
             }
             */
        }
        else
        {
            /* problems!!!! */
            printf("Error with reserving space to hold old and new photons\n");
            exit(0);
        }
        
        //also expand memory of other arrays
        tmp_double=realloc(*all_time_steps, ((*num_ph)+ph_tot-null_ph_count)*sizeof(double));
        if (tmp_double!=NULL)
        {
            *all_time_steps=tmp_double;
        }
        else
        {
            printf("Error with reallocating space to hold data about each photon's time step until an interaction occurs\n");
        }
        tmp_int=realloc(*sorted_indexes, ((*num_ph)+ph_tot-null_ph_count)*sizeof(int));
        if (tmp_int!=NULL)
        {
            *sorted_indexes=tmp_int;
        }
        else
        {
            printf("Error with reallocating space to hold data about the order in which each photon would have an interaction\n");
        }
        
        net_ph=(ph_tot-null_ph_count);
        null_ph_count=ph_tot; // use this to set the photons recently allocated as null phtoons (this can help if we decide to directly double (or *1.5) number of photons each time we need to allocate more memory, then use factor*((*num_ph)+ph_tot)-(*num_ph)
        null_ph_indexes=malloc((ph_tot+null_ph_count)*sizeof(int));
        j=0;
        for (i=((*num_ph)+net_ph)-1;i >=0 ;i--)
        {
            //fprintf(fPtr, "idx %d\n", i);
            //fflush(fPtr);
            if (((*ph_orig)[i].weight == 0)  || (i >= *num_ph))
            {
                //preset values for the the newly created spots to hold the emitted phtoons in
                (*ph_orig)[i].weight=0;
                (*ph_orig)[i].nearest_block_index=-1;
                *(null_ph_indexes+j)=i; //save this information so we can use the same syntax for both cases in saving the emitted photon data
                fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
                fflush(fPtr);
                j++;
            }
            /*
            else
            {
                //for one fo the original photon see if there are any photons with non zero weights and the nearest_block_index==-1, change the nearest_block_index to 0
                //this was an synch photon emitted in the last frame and has to have this value changes so it can be scattered/absorbed
                if (((*ph_orig)[i].weight != 0) && ((*ph_orig)[i].nearest_block_index == -1))
                {
                    //if the photon weight isnt 0 and the nearest_block_index==-1 change the nearest_block_index to 0
                    //this was an synch photon emitted in the last frame and has to have this value changes so it can be scattered/absorbed
                    (*ph_orig)[i].nearest_block_index == 0;
                    fprintf(fPtr, "Allowing photon %d to scatter/absorb now.\n", i);
                }
            }
            */
        }
        count_null_indexes=ph_tot; //use this to count the number fo null photons we have actually created, (this can help if we decide to directly double (or *1.5) number of photons each time we need to allocate more memory, then use factor*((*num_ph)+ph_tot)-(*num_ph)
        
        //loop through the original set of photons to see if
        
        fprintf(fPtr,"Val %d\n", (*(null_ph_indexes+count_null_indexes-1)));
        *num_ph+=net_ph; //update number of photons
        *num_null_ph=ph_tot-null_ph_count; //((*num_ph)+ph_tot)-(*num_ph)-ph_tot; //reserved space - emitted photons-original photons
        fprintf(fPtr,"New Num PH %d\nNew null hum_ph %d\n", *num_ph, *num_null_ph);
        fflush(fPtr);
    }
    else
    {
        //otherwise need to find the indexes of these null photons to save the newly emitted photons in them, start searching from the end of the array to efficiently find them
        //dont need to update the number of photons here
        null_ph_indexes=malloc(null_ph_count*sizeof(int));
        j=0;
        for (i=(*num_ph)-1;i>=0;i--)
        {
            /*
            if (((*ph_orig)[i].weight != 0) && ((*ph_orig)[i].nearest_block_index == -1))
            {
                //if the photon weight isnt 0 and the nearest_block_index==-1 change the nearest_block_index to 0
                //this was an synch photon emitted in the last frame and has to have this value changes so it can be scattered/absorbed
                (*ph_orig)[i].nearest_block_index == 0;
                fprintf(fPtr, "Allowing photon %d to scatter/absorb now.\n", i);
            }
            else
                */
            if ((*ph_orig)[i].weight == 0)
            {
                // if the weight is 0, this is a photons that has been absorbed and is now null
                *(null_ph_indexes+j)=i;
                j++;
                fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
                fflush(fPtr);
                
                if (j == null_ph_count)
                {
                    i=-1; //have found al the indexes and can exit the loop, dont want to do this so we can do the first part of the if statement
                }
            }
            
        }
        
        count_null_indexes=null_ph_count;
        
        *num_null_ph=null_ph_count-ph_tot;
        
    }
    
    if (inject_single_switch == 0)
    {
        //go through blocks and assign random energies/locations to proper number of photons
        ph_tot=0;
        for (i=0;i<array_length;i++)
        {
            if ((*(r+i) >= rmin)  &&   (*(r+i)  < rmax  ) && (*(theta+i)< theta_max) && (*(theta+i) >= theta_min) )
            {
                
                el_dens= (*(dens+i))/M_P;
                nu_c=calcCyclotronFreq(calcB(el_dens,*(temp+i) , epsilon_b));
                dimlesstheta=calcDimlessTheta( *(temp+i));
                max_jnu=2*jnu(nu_c/10, nu_c, dimlesstheta, el_dens);
                
                for(j=0;j<( *(ph_dens+k) ); j++ )
                {
                    //printf("flash_array_idx: %d Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e\n",i, *(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta);

                    //have to get random frequency for the photon comoving frequency
                    /*
                    y_dum=1; //initalize loop
                    yfr_dum=0;
                    while (y_dum>yfr_dum)
                    {
                        fr_dum=gsl_rng_uniform_pos(rand)*(nu_c*1e2) ;//pow(10, gsl_rng_uniform_pos(rand)*log10(nu_c*1e2)); //in Hz
                        //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
                        y_dum=gsl_rng_uniform_pos(rand)*max_jnu;
                        //printf("%lf ",fr_dum);
                        yfr_dum=jnu(fr_dum, nu_c, dimlesstheta, el_dens);
                     
                        //printf("%lf,%lf,%e \n",fr_dum, y_dum, yfr_dum);
                     
                    }
                     */
                    fr_dum=nu_c; //set the frequency directly to the cyclotron frequency
                    fprintf(fPtr, "%lf\n ",fr_dum);
                    //exit(0);
                    position_phi=gsl_rng_uniform(rand)*2*M_PI;
                    com_v_phi=gsl_rng_uniform(rand)*2*M_PI;
                    com_v_theta=gsl_rng_uniform(rand)*M_PI; //  acos((gsl_rng_uniform(rand)*2)-1) this was for compton scatt, should be isotropic now?
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
                    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
                    printf("Assigning values to struct\n");
                    
                    idx=(*(null_ph_indexes+count_null_indexes-1));
                    fprintf(fPtr, "Placing photon in index %d\n", idx);
                    (*ph_orig)[idx].p0=(*(l_boost+0));
                    (*ph_orig)[idx].p1=(*(l_boost+1));
                    (*ph_orig)[idx].p2=(*(l_boost+2));
                    (*ph_orig)[idx].p3=(*(l_boost+3));
                    (*ph_orig)[idx].comv_p0=(*(p_comv+0));
                    (*ph_orig)[idx].comv_p1=(*(p_comv+1));
                    (*ph_orig)[idx].comv_p2=(*(p_comv+2));
                    (*ph_orig)[idx].comv_p3=(*(p_comv+3));
                    //position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
                    (*ph_orig)[idx].r0= (*(x+i))*cos(position_phi); //put photons @center of the box with random phi
                    (*ph_orig)[idx].r1=(*(x+i))*sin(position_phi) ;
                    //position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
                    (*ph_orig)[idx].r2=(*(y+i)); //y coordinate in flash becomes z coordinate in MCRaT
                    (*ph_orig)[idx].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                    (*ph_orig)[idx].s1=0;
                    (*ph_orig)[idx].s2=0;
                    (*ph_orig)[idx].s3=0;
                    (*ph_orig)[idx].num_scatt=0;
                    (*ph_orig)[idx].weight=ph_weight_adjusted;
                    (*ph_orig)[idx].nearest_block_index=0; //these photons can be scattered
                    (*ph_orig)[idx].type='s';
                    //printf("%d\n",ph_tot);
                    ph_tot++; //count how many photons have been emitted
                    count_null_indexes--; //keep track fo the null photon indexes
                    
                    if ((count_null_indexes == 0) || (ph_tot == null_ph_count))
                    {
                        //if count_null_indexes is 0, then all the null photon spaces are filled with emitted photons
                        //if ph_tot is equal to what it used to be
                        i=array_length;
                        printf("Exiting Emitting loop\n");
                    }
                }
                k++;
            }
        }
    }
    else
    {
        //need to replace the scattered synch photon with another.
        //place new photon near the old one and make sure that it has the same nu_c as the other unscattered synch photons
        idx=(*(null_ph_indexes+count_null_indexes-1));
        i=(*ph_orig)[scatt_ph_index].nearest_block_index;
        
        el_dens= (*(dens+i))/M_P;
        nu_c=calcCyclotronFreq(calcB(el_dens,*(temp+i) , epsilon_b));
        
        fr_dum=nu_c; //_scatt; //set the frequency directly to the cyclotron frequency
        fprintf(fPtr, "%lf %d\n ",fr_dum, (*ph_orig)[scatt_ph_index].nearest_block_index);
        //exit(0);
        position_phi=gsl_rng_uniform(rand)*2*M_PI;
        com_v_phi=gsl_rng_uniform(rand)*2*M_PI;
        com_v_theta=gsl_rng_uniform(rand)*M_PI; //  acos((gsl_rng_uniform(rand)*2)-1) this was for compton scatt, should be isotropic now?
        
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
        lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
        
        fprintf(fPtr, "Placing photon in index %d\n", idx);
        (*ph_orig)[idx].p0=(*(l_boost+0));
        (*ph_orig)[idx].p1=(*(l_boost+1));
        (*ph_orig)[idx].p2=(*(l_boost+2));
        (*ph_orig)[idx].p3=(*(l_boost+3));
        (*ph_orig)[idx].comv_p0=(*(p_comv+0));
        (*ph_orig)[idx].comv_p1=(*(p_comv+1));
        (*ph_orig)[idx].comv_p2=(*(p_comv+2));
        (*ph_orig)[idx].comv_p3=(*(p_comv+3));
        //position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
        (*ph_orig)[idx].r0= (*(x+i))*cos(position_phi); //put photons at center of the box with random phi
        (*ph_orig)[idx].r1=(*(x+i))*sin(position_phi) ;
        //position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
        (*ph_orig)[idx].r2=(*(y+i)); //y coordinate in flash becomes z coordinate in MCRaT
        (*ph_orig)[idx].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
        (*ph_orig)[idx].s1=0;
        (*ph_orig)[idx].s2=0;
        (*ph_orig)[idx].s3=0;
        (*ph_orig)[idx].num_scatt=0;
        (*ph_orig)[idx].weight=(*ph_orig)[scatt_ph_index].weight;
        (*ph_orig)[idx].nearest_block_index=i; //these photons can be scattered
        (*ph_orig)[idx].type='s';
        
        //change position of scattered synchrotron photon to be random in the hydro grid
        position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0; //choose between -size/2 to size/2
        (*ph_orig)[scatt_ph_index].r0=(*(x+i)+position_rand)*cos(position_phi);
        (*ph_orig)[scatt_ph_index].r1=(*(x+i)+position_rand)*sin(position_phi);
        position_rand=gsl_rng_uniform_pos(rand)*(*(szx+i))-(*(szx+i))/2.0;
        (*ph_orig)[scatt_ph_index].r2=(*(y+i)+position_rand);
        
    }
    
    
    printf("(*ph_orig)[0].p0 %e (*ph_orig)[71].p0 %e (*ph_orig)[72].p0 %e (*ph_orig)[73].p0 %e\n", (*ph_orig)[0].p0, (*ph_orig)[71].p0, (*ph_orig)[96].p0, (*ph_orig)[97].p0);
    printf("At End of function\n");
    
    
    if (null_ph_count > ph_tot)
    {
        free(null_ph_indexes);
    }
    
    //exit(0);
    free(ph_dens); free(p_comv); free(boost); free(l_boost);
    //free(ph_emit);
    
    gsl_integration_workspace_free (w);
    
    return ph_tot;
    
}

int phAbsSynch(struct photon **ph_orig, int *num_ph, int *num_abs_ph, double epsilon_b, double *temp, double *dens, FILE *fPtr)
{
    //still need to deal with below issue
    //frame 213 where the absorption doesnt occur for all emitted photons and have some absorbed before/after unabsorbed photons, how to deal with this?
    //ph 97, neg lab nu in frame 210, from -1 * c/h
    int i=0, count=0, abs_ph_count=0, num_thread=omp_get_num_threads();
    double el_dens=0, nu_c=0;
    //struct photon tmp_ph;//hold temporay photon to move its data
    
    #pragma omp parallel for num_threads(num_thread) firstprivate(el_dens, nu_c) reduction(+:abs_ph_count)
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && ((*ph_orig)[i].nearest_block_index != -1))
        {
            // if the photon isnt a null photon already, see if it should be absorbed
            
            el_dens= (*(dens+(*ph_orig)[i].nearest_block_index))/M_P;
            nu_c=calcCyclotronFreq(calcB(el_dens,*(temp+(*ph_orig)[i].nearest_block_index) , epsilon_b));
            printf("photon %d has lab nu %e comv frequency %e and nu_c %e with FLASH grid number %d\n", i, (*ph_orig)[i].p0*C_LIGHT/PL_CONST, (*ph_orig)[i].comv_p0*C_LIGHT/PL_CONST, nu_c, (*ph_orig)[i].nearest_block_index);
            if (((*ph_orig)[i].comv_p0*C_LIGHT/PL_CONST <= nu_c) || ((*ph_orig)[i].type == 's'))
            {
                //if the photon has a frequency less that nu_c, it should be absorbed and becomes a null photon
                //preset values for the the newly created spots to hold the emitted phtoons in;
                
                //if this is a synchrotron photons or photons that have been scattered that were once synch photons in this frame
                printf("photon %d being absorbed\n", i);
                if (((*ph_orig)[i].type != 'i') )
                {
                    (*ph_orig)[i].weight=0;
                    (*ph_orig)[i].nearest_block_index=-1;
                    abs_ph_count+=1;
                }
                else
                {
                    //have an injected photon or previous 'c' photon that has a nu that can be absorbed
                    (*ph_orig)[i].p0=-1; //set its energy negative so we know for later analysis that it can't be used and its been absorbed, this makes it still get saves in the hdf5 files
                    (*ph_orig)[i].nearest_block_index=-1;
                    
                    //replace the potantial null photon with this photon's data
                    (*ph_orig)[count].p0=(*ph_orig)[i].p0;
                    (*ph_orig)[count].p1=(*ph_orig)[i].p1;
                    (*ph_orig)[count].p2=(*ph_orig)[i].p2;
                    (*ph_orig)[count].p3=(*ph_orig)[i].p3;
                    (*ph_orig)[count].comv_p0=(*ph_orig)[i].comv_p0;
                    (*ph_orig)[count].comv_p1=(*ph_orig)[i].comv_p1;
                    (*ph_orig)[count].comv_p2=(*ph_orig)[i].comv_p2;
                    (*ph_orig)[count].comv_p3=(*ph_orig)[i].comv_p3;
                    (*ph_orig)[count].r0= (*ph_orig)[i].r0;
                    (*ph_orig)[count].r1=(*ph_orig)[i].r1 ;
                    (*ph_orig)[count].r2=(*ph_orig)[i].r2;
                    (*ph_orig)[count].s0=(*ph_orig)[i].s0;
                    (*ph_orig)[count].s1=(*ph_orig)[i].s1;
                    (*ph_orig)[count].s2=(*ph_orig)[i].s2;
                    (*ph_orig)[count].s3=(*ph_orig)[i].s3;
                    (*ph_orig)[count].num_scatt=(*ph_orig)[i].num_scatt;
                    (*ph_orig)[count].weight=(*ph_orig)[i].weight;
                    (*ph_orig)[count].nearest_block_index=(*ph_orig)[i].nearest_block_index;
                    (*ph_orig)[count].type=(*ph_orig)[i].type;
                    
                    count+=1; //increment count (counts non-null photons in array)
                    
                }
            }
            else
            {
                //if the phootn isnt going to be absorbed, see if its a 'c' photon thats survived and change it to an injected type
                if (((*ph_orig)[i].type == 'c') )
                {
                    (*ph_orig)[i].type = 'i';
                }
                
                //replace the potantial null photon with this photon's data
                (*ph_orig)[count].p0=(*ph_orig)[i].p0;
                (*ph_orig)[count].p1=(*ph_orig)[i].p1;
                (*ph_orig)[count].p2=(*ph_orig)[i].p2;
                (*ph_orig)[count].p3=(*ph_orig)[i].p3;
                (*ph_orig)[count].comv_p0=(*ph_orig)[i].comv_p0;
                (*ph_orig)[count].comv_p1=(*ph_orig)[i].comv_p1;
                (*ph_orig)[count].comv_p2=(*ph_orig)[i].comv_p2;
                (*ph_orig)[count].comv_p3=(*ph_orig)[i].comv_p3;
                (*ph_orig)[count].r0= (*ph_orig)[i].r0;
                (*ph_orig)[count].r1=(*ph_orig)[i].r1 ;
                (*ph_orig)[count].r2=(*ph_orig)[i].r2;
                (*ph_orig)[count].s0=(*ph_orig)[i].s0;
                (*ph_orig)[count].s1=(*ph_orig)[i].s1;
                (*ph_orig)[count].s2=(*ph_orig)[i].s2;
                (*ph_orig)[count].s3=(*ph_orig)[i].s3;
                (*ph_orig)[count].num_scatt=(*ph_orig)[i].num_scatt;
                (*ph_orig)[count].weight=(*ph_orig)[i].weight;
                (*ph_orig)[count].nearest_block_index=(*ph_orig)[i].nearest_block_index;
                (*ph_orig)[count].type=(*ph_orig)[i].type;
                
                //increment count
                count+=1;
                
                
            }
        }
        else
        {
            //see if the photon was a previous 'i' photon absorbed that we still have to account for in the array
            if (((*ph_orig)[i].p0 < 0) )
            {
                //replace the potantial null photon with this photon's data
                (*ph_orig)[count].p0=(*ph_orig)[i].p0;
                (*ph_orig)[count].p1=(*ph_orig)[i].p1;
                (*ph_orig)[count].p2=(*ph_orig)[i].p2;
                (*ph_orig)[count].p3=(*ph_orig)[i].p3;
                (*ph_orig)[count].comv_p0=(*ph_orig)[i].comv_p0;
                (*ph_orig)[count].comv_p1=(*ph_orig)[i].comv_p1;
                (*ph_orig)[count].comv_p2=(*ph_orig)[i].comv_p2;
                (*ph_orig)[count].comv_p3=(*ph_orig)[i].comv_p3;
                (*ph_orig)[count].r0= (*ph_orig)[i].r0;
                (*ph_orig)[count].r1=(*ph_orig)[i].r1 ;
                (*ph_orig)[count].r2=(*ph_orig)[i].r2;
                (*ph_orig)[count].s0=(*ph_orig)[i].s0;
                (*ph_orig)[count].s1=(*ph_orig)[i].s1;
                (*ph_orig)[count].s2=(*ph_orig)[i].s2;
                (*ph_orig)[count].s3=(*ph_orig)[i].s3;
                (*ph_orig)[count].num_scatt=(*ph_orig)[i].num_scatt;
                (*ph_orig)[count].weight=(*ph_orig)[i].weight;
                (*ph_orig)[count].nearest_block_index=(*ph_orig)[i].nearest_block_index;
                (*ph_orig)[count].type=(*ph_orig)[i].type;
                
                //increment count
                count+=1;
            }
        }
        
        fprintf(fPtr, "photon %d has lab frequency %e and weight %e with FLASH grid number %d\n", i, (*ph_orig)[i].p0*C_LIGHT/PL_CONST, (*ph_orig)[i].weight, (*ph_orig)[i].nearest_block_index);
    }
    *num_abs_ph=abs_ph_count;
    /*
    count=0;
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && ((*ph_orig)[i].nearest_block_index != -1))
        {
            fprintf(fPtr, "photon %d has lab frequency %e and weight %e with FLASH grid number %d\n", i, (*ph_orig)[i].p0*C_LIGHT/PL_CONST, (*ph_orig)[i].weight, (*ph_orig)[i].nearest_block_index);
            
            tmp_ph.p0=(*ph_orig)[i].p0;
            tmp_ph.p1=(*ph_orig)[i].p1;
            tmp_ph.p2=(*ph_orig)[i].p2;
            tmp_ph.p3=(*ph_orig)[i].p3;
            tmp_ph.comv_p0=(*ph_orig)[i].comv_p0;
            tmp_ph.comv_p1=(*ph_orig)[i].comv_p1;
            tmp_ph.comv_p2=(*ph_orig)[i].comv_p2;
            tmp_ph.comv_p3=(*ph_orig)[i].comv_p3;
            tmp_ph.r0= (*ph_orig)[i].r0;
            tmp_ph.r1=(*ph_orig)[i].r1 ;
            tmp_ph.r2=(*ph_orig)[i].r2;
            tmp_ph.s0=(*ph_orig)[i].s0;
            tmp_ph.s1=(*ph_orig)[i].s1;
            tmp_ph.s2=(*ph_orig)[i].s2;
            tmp_ph.s3=(*ph_orig)[i].s3;
            tmp_ph.num_scatt=(*ph_orig)[i].num_scatt;
            tmp_ph.weight=(*ph_orig)[i].weight;
            tmp_ph.nearest_block_index=(*ph_orig)[i].nearest_block_index;
            tmp_ph.type=(*ph_orig)[i].type;
            
            //replace the potantial null photon with this photon's data
            (*ph_orig)[count].p0=tmp_ph.p0;
            (*ph_orig)[count].p1=tmp_ph.p1;
            (*ph_orig)[count].p2=tmp_ph.p2;
            (*ph_orig)[count].p3=tmp_ph.p3;
            (*ph_orig)[count].comv_p0=tmp_ph.comv_p0;
            (*ph_orig)[count].comv_p1=tmp_ph.comv_p1;
            (*ph_orig)[count].comv_p2=tmp_ph.comv_p2;
            (*ph_orig)[count].comv_p3=tmp_ph.comv_p3;
            (*ph_orig)[count].r0= tmp_ph.r0;
            (*ph_orig)[count].r1=tmp_ph.r1 ;
            (*ph_orig)[count].r2=tmp_ph.r2;
            (*ph_orig)[count].s0=tmp_ph.s0;
            (*ph_orig)[count].s1=tmp_ph.s1;
            (*ph_orig)[count].s2=tmp_ph.s2;
            (*ph_orig)[count].s3=tmp_ph.s3;
            (*ph_orig)[count].num_scatt=tmp_ph.num_scatt;
            (*ph_orig)[count].weight=tmp_ph.weight;
            (*ph_orig)[count].nearest_block_index=tmp_ph.nearest_block_index;
            (*ph_orig)[count].type=tmp_ph.type;
            
            count+=1; //increment count (counts non-null photons in array)
        }
    }
    */
    while (count<*num_ph)
    {
        //overwrite the last few photons to make sure that they are null photons
        (*ph_orig)[count].weight=0;
        (*ph_orig)[count].nearest_block_index=-1;
        fprintf(fPtr, "photon %d has frequency %e and weight %e with FLASH grid number %d\n", count, (*ph_orig)[count].comv_p0*C_LIGHT/PL_CONST, (*ph_orig)[count].weight, (*ph_orig)[count].nearest_block_index);
        fflush(fPtr);
        
        count+=1;
    }
    
    return 0;
}
