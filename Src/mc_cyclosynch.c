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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_integration.h>

double calcCyclotronFreq(double magnetic_field)
{
    //B has to be in gauss
    return CHARGE_EL*magnetic_field/(2*M_PI*M_EL*C_LIGHT);
}

double calcEB(double magnetic_field)
{
    //helper function to compare to vurm 2013
    return PL_CONST*calcCyclotronFreq(magnetic_field);
}

double calcBoundaryE(double magnetic_field, double temp)
{
    //helper function to compare to vurm 2013
    return 14*pow(M_EL*C_LIGHT*C_LIGHT, 1.0/10.0)*pow(calcEB(magnetic_field), 9.0/10.0)*pow(calcDimlessTheta(temp), 3.0/10.0);
}

double calcDimlessTheta(double temp)
{
    //temp has to be in kelvin
    return K_B*temp/(M_EL*C_LIGHT*C_LIGHT);
}

double calcB(double el_dens, double temp)
{
    //calc the B field from assuming its some fraction of the matter energy density
    //assume equipartition here
    #if B_FIELD_CALC == INTERNAL_E
        return sqrt(EPSILON_B*8*M_PI*3*el_dens*K_B*temp/2);
    #else
        //otherwise calculate B from the total energy
        return sqrt(8*M_PI*EPSILON_B*(el_dens*M_P*C_LIGHT*C_LIGHT+4*A_RAD*temp*temp*temp*temp/3));
    #endif
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
    
    return (8*M_PI*nu*nu)/(exp(PL_CONST*nu/(K_B*temp))-1)/(C_LIGHT*C_LIGHT*C_LIGHT);
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

double synCrossSection(double el_dens, double T, double nu_ph, double p_el)
{
    double b_cr=FINE_STRUCT*sqrt(M_EL*C_LIGHT*C_LIGHT/pow(R_EL,3.0));
    double B=calcB(el_dens, T);
    double nu_c=calcCyclotronFreq(B);
    double gamma_el=sqrt(p_el*p_el+1);
    
    //printf("calc gamma %e, temp %e\n", gamma_el, T);
    
    return (3.0*M_PI*M_PI/8.0)*(THOM_X_SECT/FINE_STRUCT)*(b_cr/B)*pow(nu_c/nu_ph, 2.0) * exp(-2*nu_ph*(gamma_el*log((gamma_el+1)/p_el)-1)/nu_c)* ((C(nu_ph, nu_c, gamma_el, p_el)/G(gamma_el, p_el))-(G_prime(gamma_el, p_el)/pow(G(gamma_el, p_el),2.0)));
}

double calcCyclosynchRLimits(int frame_scatt, int frame_inj, double fps,  double r_inj, char *min_or_max)
{
    double val=r_inj;
    if (strcmp(min_or_max, "min")==0)
    {
        //printf("IN MIN\nframe_scatt %e frame_inj %e fps %e r_inj %e C_LIGHT %e\n", frame_scatt, frame_inj, fps, r_inj, C_LIGHT);
        val+=(C_LIGHT*(frame_scatt-frame_inj)/fps - 0.5*C_LIGHT/fps);
    }
    else
    {
        //printf("IN MAX\n");
        val+=(C_LIGHT*(frame_scatt-frame_inj)/fps + 0.5*C_LIGHT/fps);
    }
    
    //printf("Val %e\n", val);
    
    return val;
}

int rebinCyclosynchCompPhotons(struct photon **ph_orig, int *num_ph,  int *num_null_ph, int *num_cyclosynch_ph_emit, int *scatt_cyclosynch_num_ph, double **all_time_steps, int **sorted_indexes, int max_photons, double thread_theta_min, double thread_theta_max , gsl_rng * rand, FILE *fPtr)
{
    int i=0, j=0, k=0, count=0, count_c_ph=0, end_count=(*scatt_cyclosynch_num_ph), idx=0, num_thread=1;
    #if defined(_OPENMP)
    num_thread=omp_get_num_threads();
    #endif
    int synch_comp_photon_count=0, synch_photon_count=0, num_avg=12, num_bins=(0.1)*max_photons; //some factor of the max number of photons that is specified in the mc.par file, num bins is also in test function
    double dtheta_bin=0.5*M_PI/180; //the size of the bin that we want to produce for spatial binning in theta
    int num_bins_theta=(thread_theta_max-thread_theta_min)/dtheta_bin;//try this many bins such that we have 0.5 degree resolution, can also try to do adaptive binning with constant SNR
    double avg_values[12]={0}; //number of averages that'll be taken is given by num_avg in above line
    double p0_min=DBL_MAX, p0_max=0, log_p0_min=0, log_p0_max=0;//look at p0 of photons not by frequency since its just nu=p0*C_LIGHT/PL_CONST
    double rand1=0, rand2=0, phi=0, theta=0;
    double min_range=0, max_range=0, energy=0;
    double ph_r=0, ph_theta=0, temp_theta_max=0, temp_theta_min=DBL_MAX;
    //int *synch_comp_photon_idx=NULL; make this an array b/c had issue with deallocating this memory for some reason
    int synch_comp_photon_idx[*scatt_cyclosynch_num_ph];
    //struct photon *rebin_ph=malloc(num_bins* sizeof (struct photon ));
    int num_null_rebin_ph=0;
    struct photon *tmp=NULL;
    double *tmp_double=NULL;
    int *tmp_int=NULL;
    gsl_histogram * h = gsl_histogram_alloc (num_bins);
    gsl_histogram2d * h_phi_theta = gsl_histogram2d_alloc (360, 180); //x is for phi goes from 0 to 2pi and y is for theta, goes from 0 to pi
    gsl_histogram2d_set_ranges_uniform (h_phi_theta, 0.0, 360.0,0.0, 180.0); //set the ranges to be 1 degree wide
    gsl_histogram2d_pdf * pdf_phi_theta = gsl_histogram2d_pdf_alloc (h_phi_theta->nx, h_phi_theta->ny);
    
    //calc min and max p0 and set the bins to be even within this interval
    //save the frequencies of photons to a
    //#pragma omp parallel for num_threads(num_thread) reduction(min:nu_min) reduction(max:nu_max)
    //synch_comp_photon_idx=malloc((*scatt_cyclosynch_num_ph)*sizeof(int));
    //if (synch_comp_photon_idx==NULL)
    //{
    //    printf("Error with allocating space to hold data for synch_comp_photon_idx\n");
    //    exit(1);
    //}
    
    
    fprintf(fPtr, "In the rebin func; num_threads %d scatt_cyclosynch_num_ph %d, num_ph %d\n", num_thread, (*scatt_cyclosynch_num_ph), *num_ph);
    fflush(fPtr); //try to see if all the photons are being set to some value at the end of this function when allocating things
    
    /*
    count=0;
    for (i=0;i<*num_ph;i++)
    {
        //fprintf(fPtr, "%d %c %e %e\n", i, (*ph_orig)[i].type, (*ph_orig)[i].weight, (*ph_orig)[i].p0 );
        //fflush(fPtr);
        
        if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
        {
            count++;
        }
    }
    fprintf(fPtr, "in rebin: count is: %d and scatt_cyclosynch_num_ph is %d\n", count, *scatt_cyclosynch_num_ph );
    */
    int min_idx=0, max_idx=0;
    count=0;
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
        {
            //see if the photon's nu is larger than nu_max or smaller than nu_min
            if (((*ph_orig)[i].p0< p0_min))
            {
                //dont include any absorbed UNABSORBED_CS_PHOTON photons that have negative P0 values
                p0_min= (*ph_orig)[i].p0;
                min_idx=i;
                //fprintf(fPtr, "new p0 min %e\n", (p0_min) );
            }
            
            if ((*ph_orig)[i].p0> p0_max)
            {
                p0_max= (*ph_orig)[i].p0;
                max_idx=i;
                //fprintf(fPtr, "new p0 max %e\n", (p0_max) );
            }
            
            //look at min and max theta of photons
            ph_r=pow(((*ph_orig)[i].r0)*((*ph_orig)[i].r0) + ((*ph_orig)[i].r1)*((*ph_orig)[i].r1) + ((*ph_orig)[i].r2)*((*ph_orig)[i].r2),0.5);
            ph_theta=acos(((*ph_orig)[i].r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            
            if (ph_theta > temp_theta_max )
            {
                temp_theta_max=ph_theta;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_theta<temp_theta_min)
            {
                temp_theta_min=ph_theta;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }

            
            // also save the index of these photons because they wil become null later on
            //*(synch_comp_photon_idx+count)=i;
            synch_comp_photon_idx[count]=i;
            //fprintf(fPtr, "Save index %d\n", i );
            count++;
            
            if ((*ph_orig)[i].type == COMPTONIZED_PHOTON)
            {
                //keep track of the number of COMPTONIZED_PHOTON photons so we can know if the array needs to be increased in size, also take num_null_ph into account in doing this
                count_c_ph+=1;
            }
        }
        else if (((*ph_orig)[i].type == CS_POOL_PHOTON) && ((*ph_orig)[i].weight != 0))
        {
            synch_photon_count++;
        }
    }
    
    
    //temp_theta_min=floor(temp_theta_min*180/M_PI)*M_PI/180;
    //temp_theta_max=ceil(temp_theta_max*180/M_PI)*M_PI/180;
    num_bins_theta=1+(temp_theta_max-temp_theta_min)/dtheta_bin;
    
    fprintf(fPtr, "min, max (keV): %e %e log p0 min, max: %e %e idx: %d %d\n", p0_min*C_LIGHT/1.6e-9,p0_max*C_LIGHT/1.6e-9 , log10(p0_min), log10(p0_max), min_idx, max_idx );
    fprintf(fPtr, "min, max (theta in deg): %e %e number of bins %d count: %d\n", temp_theta_min*180/M_PI, temp_theta_max*180/M_PI, num_bins_theta, count );
    fflush(fPtr);
    
    if (num_bins_theta*num_bins>=max_photons)
    {
        fprintf(fPtr, "The number of rebinned photons, %d, is larger than max_photons %d nd will not rebin efficiently. Adjust the parameters such that the number of bins in theta and energy are less than the number of photons that will lead to rebinning.\n",  num_bins_theta*num_bins, max_photons);
        fflush(fPtr);
        printf( "The number of rebinned photons, %d, is larger than max_photons %d nd will not rebin efficiently. Adjust the parameters such that the number of bins in theta and energy are less than the number of photons that will lead to rebinning.\n",  num_bins_theta*num_bins, max_photons);
        exit(1);

    }

    gsl_histogram_set_ranges_uniform (h, log10(p0_min), log10(p0_max*(1+1e-6)));
    
    gsl_histogram * h_theta_space = gsl_histogram_alloc (num_bins_theta);
    gsl_histogram_set_ranges_uniform (h_theta_space, temp_theta_min, temp_theta_max+dtheta_bin);
    
    struct photon *rebin_ph=malloc(num_bins*num_bins_theta* sizeof (struct photon ));
    
    gsl_histogram2d * h_energy_theta = gsl_histogram2d_alloc (num_bins, num_bins_theta); //x is for energy  and y is for spatial theta, goes from 0 to pi
    gsl_histogram2d_set_ranges_uniform (h_energy_theta, log10(p0_min), log10(p0_max*(1+1e-6)), temp_theta_min, temp_theta_max+dtheta_bin);



    
    //populate histogram for photons with nu that falss within the proper histogram bin
    //may not need this loop, can just check if the photon nu falls within the bin edges and do averages etc within next loop
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
        {
            //gsl_histogram_accumulate (h, log10((*ph_orig)[i].p0), (*ph_orig)[i].weight);
            gsl_histogram_increment (h, log10((*ph_orig)[i].p0));
            
            ph_r=pow(((*ph_orig)[i].r0)*((*ph_orig)[i].r0) + ((*ph_orig)[i].r1)*((*ph_orig)[i].r1) + ((*ph_orig)[i].r2)*((*ph_orig)[i].r2),0.5);
            ph_theta=acos(((*ph_orig)[i].r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            gsl_histogram_increment (h_theta_space, ph_theta);
            
            gsl_histogram2d_increment(h_energy_theta, log10((*ph_orig)[i].p0), ph_theta);

        }
    }
    
    //gsl_histogram_fprintf(fPtr, h_theta_space, "%g", "%g");
    gsl_histogram2d_fprintf(fPtr, h_energy_theta, "%g", "%g");

    //for the photons that fall within a given nu bin, histogram thier theta and phi and choose a random number to sample from the distribution to get the new photons' 4 momentum, to parallelize this can put num+ph lop outside and cpunt loop inside and make avg_value array 2D with count index and other index being the average values
    
    for (count=0;count<num_bins;count++)
    {

        if (gsl_histogram_get(h, count) > 1)
        {
            for (j=0;j<num_avg;j++)
            {
                avg_values[j]=0;
            }
            
            //loop over the number of photons
            for (i=0;i<*num_ph;i++)
            {
                if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
                {
                    gsl_histogram_get_range(h, count, &min_range, &max_range);
                    //if the photon nu falls in the count bin of the nu histogram then add it to the phi_theta 2d hist
                    if ((log10((*ph_orig)[i].p0)< max_range  ) && (log10((*ph_orig)[i].p0)>=min_range))
                    {
                        gsl_histogram2d_increment(h_phi_theta, fmod(atan2((*ph_orig)[i].p2,((*ph_orig)[i].p1)*180/M_PI + 360),360.0), (180/M_PI)*acos(((*ph_orig)[i].p3)/((*ph_orig)[i].p0)) );
                        
                        avg_values[0] += (*ph_orig)[i].r0*(*ph_orig)[i].weight; //used to calc weighted averages
                        avg_values[1] += (*ph_orig)[i].r1*(*ph_orig)[i].weight;
                        avg_values[2] += (*ph_orig)[i].r2*(*ph_orig)[i].weight;
                        avg_values[3] += (*ph_orig)[i].s0*(*ph_orig)[i].weight;
                        avg_values[4] += (*ph_orig)[i].s1*(*ph_orig)[i].weight;
                        avg_values[5] += (*ph_orig)[i].s2*(*ph_orig)[i].weight;
                        avg_values[6] += (*ph_orig)[i].s3*(*ph_orig)[i].weight;
                        avg_values[7] += (*ph_orig)[i].num_scatt*(*ph_orig)[i].weight;
                        avg_values[8] += (*ph_orig)[i].weight;
                        
                        //average theta and phi of photons
                        avg_values[9] += fmod(atan2((*ph_orig)[i].p2,((*ph_orig)[i].p1)*180/M_PI + 360),360.0) *(*ph_orig)[i].weight;
                        avg_values[10] += (180/M_PI)*acos(((*ph_orig)[i].p3)/((*ph_orig)[i].p0))*(*ph_orig)[i].weight;
                        avg_values[11] +=(*ph_orig)[i].p0*(*ph_orig)[i].weight;
                    }
                }
            }
            
            //fprintf(fPtr, "bin %e-%e has %e photons\n",pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, gsl_histogram_get(h, count));
            
            energy=avg_values[11]/avg_values[8];//pow(10,0.5*(max_range+min_range));
            
            //initiate pdf as the histogram of phi and theta
            gsl_histogram2d_pdf_init (pdf_phi_theta, h_phi_theta);
            
            //get two random values
            rand1=gsl_rng_uniform(rand);
            rand2=gsl_rng_uniform(rand);
            
            //choose random phi and theta value
            gsl_histogram2d_pdf_sample (pdf_phi_theta, rand1, rand2, &phi, &theta);//phi and theta are in degreesneed to convert into radians later
            
            phi=avg_values[9]/avg_values[8];
            theta=avg_values[10]/avg_values[8];
            
            (rebin_ph+count)->type = COMPTONIZED_PHOTON;
            
            (rebin_ph+count)->p0=energy;
            (rebin_ph+count)->p1=energy*sin(theta*M_PI/180)*cos(phi*M_PI/180);
            (rebin_ph+count)->p2=energy*sin(theta*M_PI/180)*sin(phi*M_PI/180);
            (rebin_ph+count)->p3=energy*cos(theta*M_PI/180);
            (rebin_ph+count)->comv_p0=0;
            (rebin_ph+count)->comv_p1=0;
            (rebin_ph+count)->comv_p2=0;
            (rebin_ph+count)->comv_p3=0;
            (rebin_ph+count)->r0=avg_values[0]/avg_values[8];
            (rebin_ph+count)->r1= avg_values[1]/avg_values[8];
            (rebin_ph+count)->r2=avg_values[2]/avg_values[8];
            (rebin_ph+count)->s0=avg_values[3]/avg_values[8]; // stokes parameterized are normalized such that I always =1
            (rebin_ph+count)->s1=avg_values[4]/avg_values[8];
            (rebin_ph+count)->s2=avg_values[5]/avg_values[8];
            (rebin_ph+count)->s3=avg_values[6]/avg_values[8];
            (rebin_ph+count)->num_scatt=avg_values[7]/avg_values[8];
            (rebin_ph+count)->weight=avg_values[8];
            (rebin_ph+count)->nearest_block_index=0; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
                    
            
            
            //gsl_histogram2d_fprintf (stdout, h_phi_theta, "%g", "%g");
            //fprintf(fPtr, "Chosen phi: %e chosen theta: %e weight: %e\n\n", phi, theta, avg_values[8] );
            //reset the histogram and the pdf
            gsl_histogram2d_reset(pdf_phi_theta);
            gsl_histogram_reset(h_phi_theta);
        }
        else if (gsl_histogram_get(h, count) == 1)
        {
            gsl_histogram_get_range(h, count, &min_range, &max_range);
            //fprintf(fPtr, "bin %e-%e has %e photons\n", pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, 1.0);

            //for thr case of just 1 hoton being in the bin just set the rebinned photon to the one photons parameters
            for (i=0;i<*num_ph;i++)
            {
                if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON))&& ((*ph_orig)[i].p0 > 0))
                {
                    if ((log10((*ph_orig)[i].p0)< max_range  ) && (log10((*ph_orig)[i].p0)>=min_range))
                    {
                        (rebin_ph+count)->p0=(*ph_orig)[i].p0;
                        (rebin_ph+count)->p1=(*ph_orig)[i].p1;
                        (rebin_ph+count)->p2=(*ph_orig)[i].p2;
                        (rebin_ph+count)->p3=(*ph_orig)[i].p3;
                        (rebin_ph+count)->comv_p0=(*ph_orig)[i].comv_p0;
                        (rebin_ph+count)->comv_p1=(*ph_orig)[i].comv_p1;
                        (rebin_ph+count)->comv_p2=(*ph_orig)[i].comv_p2;
                        (rebin_ph+count)->comv_p3=(*ph_orig)[i].comv_p3;
                        (rebin_ph+count)->r0=(*ph_orig)[i].r0;
                        (rebin_ph+count)->r1= (*ph_orig)[i].r1;
                        (rebin_ph+count)->r2=(*ph_orig)[i].r2; //y coordinate in flash becomes z coordinate in MCRaT
                        (rebin_ph+count)->s0=(*ph_orig)[i].s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                        (rebin_ph+count)->s1=(*ph_orig)[i].s1;
                        (rebin_ph+count)->s2=(*ph_orig)[i].s2;
                        (rebin_ph+count)->s3=(*ph_orig)[i].s3;
                        (rebin_ph+count)->num_scatt=(*ph_orig)[i].num_scatt;
                        (rebin_ph+count)->weight=(*ph_orig)[i].weight;
                        (rebin_ph+count)->nearest_block_index=(*ph_orig)[i].nearest_block_index; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
                        
                        
                        i=*num_ph;
                    }
                }
            }
            
        }
        else
        {
            //fprintf(fPtr, "Rebinned Photon is a null photon because there are no photons in this energy bin.\n");
            (rebin_ph+count)->type = COMPTONIZED_PHOTON;
            
            gsl_histogram_get_range(h, count, &min_range, &max_range);
            energy=pow(10,0.5*(max_range+min_range));
            
            //fprintf(fPtr, "bin %e-%e has %e photons\n", pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, 0.0);

            
            (rebin_ph+count)->p0=energy;
            (rebin_ph+count)->p1=0;
            (rebin_ph+count)->p2=0;
            (rebin_ph+count)->p3=0;
            (rebin_ph+count)->comv_p0=0;
            (rebin_ph+count)->comv_p1=0;
            (rebin_ph+count)->comv_p2=0;
            (rebin_ph+count)->comv_p3=0;
            (rebin_ph+count)->r0=0;
            (rebin_ph+count)->r1= 0;
            (rebin_ph+count)->r2=0;
            (rebin_ph+count)->s0=1; // stokes parameterized are normalized such that I always =1
            (rebin_ph+count)->s1=0;
            (rebin_ph+count)->s2=0;
            (rebin_ph+count)->s3=0;
            (rebin_ph+count)->num_scatt=0;
            (rebin_ph+count)->weight=0;
            (rebin_ph+count)->nearest_block_index=-1; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function

        }
        
        
    }
    
    
    //find indexes of old photons that will not become null photons
    //if the photons are UNABSORBED_CS_PHOTON photons make them have p0=-1
    //for the COMPTONIZED_PHOTON photons replace the first num_bins indexes with the new rebinned photons and replace the rest of the indexes with null values
    //may need to expand the array of photons, shouldnt need to do this though
    
    //this is a default setting, see comment below where there was an else statement that had the failing line
    //end_count=(scatt_cyclosynch_num_ph);
    
    if ((count_c_ph+(*num_null_ph))<num_bins)
    {
        //need to expand the array
        //if the totoal number of photons to be emitted is larger than the number of null phtons curently in the array, then have to grow the array
        //need to realloc memory to hold the old photon info and the new emitted photon's info
        //before was doing ((*num_ph)+num_bins-count_c_ph-(*num_null_ph) ) but instead did the beow since that creatd extr space for null photons to be filled in later on
        fprintf(fPtr, "Rebin: Allocating %d space\n", ((*num_ph)+num_bins-count_c_ph+(*num_null_ph) )); //befoe was dong
        fflush(fPtr);
        tmp=realloc(*ph_orig, ((*num_ph)+num_bins-count_c_ph+(*num_null_ph))* sizeof (struct photon )); //may have to look into directly doubling (or *1.5) number of photons each time we need to allocate more memory, can do after looking at profiling for "just enough" memory method
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
            exit(1);
        }
        //fprintf(fPtr, "Rebin: Allocating %d space good in 1st realloc\n", ((*num_ph)+num_bins-count_c_ph+(*num_null_ph) ));
        //fflush(fPtr);

        
        //also expand memory of other arrays
        tmp_double=realloc(*all_time_steps, ((*num_ph)+num_bins-count_c_ph+(*num_null_ph))*sizeof(double));
        if (tmp_double!=NULL)
        {
            *all_time_steps=tmp_double;
        }
        else
        {
            printf("Error with reallocating space to hold data about each photon's time step until an interaction occurs\n");
            exit(1);

        }
        
        //fprintf(fPtr, "Rebin: Allocating %d space good in 2nd realloc\n", ((*num_ph)+num_bins-count_c_ph+(*num_null_ph) ));
        //fflush(fPtr);

        
        tmp_int=realloc(*sorted_indexes, ((*num_ph)+num_bins-count_c_ph+(*num_null_ph))*sizeof(int));
        if (tmp_int!=NULL)
        {
            *sorted_indexes=tmp_int;
        }
        else
        {
            printf("Error with reallocating space to hold data about the order in which each photon would have an interaction\n");
            exit(1);

        }
        //fprintf(fPtr, "Rebin: Allocating %d space good in 3rd realloc\n", ((*num_ph)+num_bins-count_c_ph+(*num_null_ph) ));
        //fflush(fPtr);

        //tmp_int=realloc(*synch_comp_photon_idx, ((*scatt_cyclosynch_num_ph)+num_bins-count_c_ph-(*num_null_ph))*sizeof(int));
        //if (tmp_int!=NULL)
        //{
        //    *synch_comp_photon_idx=tmp_int;
        //}
        //else
        //{
        //    printf("Error with reallocating space to hold data about the order in which each photon would have an interaction\n");
        //}
        /*
        net_ph=num_bins-count_c_ph+(*num_null_ph);
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
         
        }
        count_null_indexes=ph_tot; //use this to count the number fo null photons we have actually created, (this can help if we decide to directly double (or *1.5) number of photons each time we need to allocate more memory, then use factor*((*num_ph)+ph_tot)-(*num_ph)
        
        //loop through the original set of photons to see if
        
        fprintf(fPtr,"Val %d\n", (*(null_ph_indexes+count_null_indexes-1)));
        *num_ph+=net_ph; //update number of photons
        *num_null_ph=ph_tot-null_ph_count; //((*num_ph)+ph_tot)-(*num_ph)-ph_tot; //reserved space - emitted photons-original photons
        fprintf(fPtr,"old Num PH %d\n", *num_ph);
        fflush(fPtr);
         */
        
        *num_ph=( *num_ph)+(num_bins-count_c_ph+(*num_null_ph));
        end_count=(*scatt_cyclosynch_num_ph)+num_bins-count_c_ph+(*num_null_ph);
    }
    //else dont knwo why this is failing try to put this before if, so it gets sets and only if the above if statement is true will end_count be modified
    //{
    //    end_count=(*scatt_cyclosynch_num_ph);
    //}
    
    //go through and assign the rebinned photons to the COMPTONIZED_PHOTON phtoons and make all UNABSORBED_CS_PHOTON photons become "absorbed"
    j=0;
    count=0;
    i=0;
    for (i=0;i<end_count;i++)
    {
        if (i<(*scatt_cyclosynch_num_ph))
        {
            //the photon idx can be found from the original array
            //idx=(*(synch_comp_photon_idx+i));
            idx=synch_comp_photon_idx[i];
        }
        else
        {
            //enter this if we realloc the arrays
            //idx=(*num_ph)+i;
            idx=i-(*scatt_cyclosynch_num_ph)+( *num_ph)-(num_bins-count_c_ph+(*num_null_ph));//go to the end of the old num_ph value and start setting things to null photons
        }
        /*
        if (idx==1183)
        {
            printf("HERE IN REBIN\n");
        }
        */
        /*
        if ((*ph_orig)[idx].type == UNABSORBED_CS_PHOTON)
        {
            //if the photon is a compton scatered synch photon from a past frame treat it as though it has been absorbed
            (*ph_orig)[idx].p0=-1; //set its energy negative so we know for later analysis that it can't be used and its been "absorbed", this makes it still get saves in the hdf5 files
            (*ph_orig)[idx].nearest_block_index=-1;
            (*ph_orig)[idx].weight=0; //now making the absorbed UNABSORBED_CS_PHOTON photons become null photons
        }
        else if ((*ph_orig)[idx].type == COMPTONIZED_PHOTON)
         */
        if (((*ph_orig)[idx].type == UNABSORBED_CS_PHOTON) || ((*ph_orig)[idx].type == COMPTONIZED_PHOTON))
        {
            if (count<num_bins)
            {
                //keep saving the rebinned photns
                (*ph_orig)[idx].p0=(rebin_ph+count)->p0;
                (*ph_orig)[idx].p1=(rebin_ph+count)->p1;
                (*ph_orig)[idx].p2=(rebin_ph+count)->p2;
                (*ph_orig)[idx].p3=(rebin_ph+count)->p3;
                (*ph_orig)[idx].comv_p0=(rebin_ph+count)->comv_p0;
                (*ph_orig)[idx].comv_p1=(rebin_ph+count)->comv_p1;
                (*ph_orig)[idx].comv_p2=(rebin_ph+count)->comv_p2;
                (*ph_orig)[idx].comv_p3=(rebin_ph+count)->comv_p3;
                (*ph_orig)[idx].r0=(rebin_ph+count)->r0;
                (*ph_orig)[idx].r1=(rebin_ph+count)->r1;
                (*ph_orig)[idx].r2=(rebin_ph+count)->r2; //y coordinate in flash becomes z coordinate in MCRaT
                (*ph_orig)[idx].s0=(rebin_ph+count)->s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                (*ph_orig)[idx].s1=(rebin_ph+count)->s1;
                (*ph_orig)[idx].s2=(rebin_ph+count)->s2;
                (*ph_orig)[idx].s3=(rebin_ph+count)->s3;
                (*ph_orig)[idx].num_scatt=(rebin_ph+count)->num_scatt;
                (*ph_orig)[idx].weight=(rebin_ph+count)->weight;
                (*ph_orig)[idx].nearest_block_index=(rebin_ph+count)->nearest_block_index;
                (*ph_orig)[idx].type = COMPTONIZED_PHOTON;
                
                
                if ((rebin_ph+count)->weight==0)
                {
                    //if the bin had no photons in it, the rebinned photon is effectively null
                    //j++;
                    num_null_rebin_ph++;
                }
                
                count++;
                
            }
            else
            {
                if ((*ph_orig)[idx].type == UNABSORBED_CS_PHOTON)
                {
                    //if the photon is a compton scatered synch photon from a past frame treat it as though it has been absorbed
                    (*ph_orig)[idx].p0=-1; //set its energy negative so we know for later analysis that it can't be used and its been "absorbed", this makes it still get saves in the hdf5 files
                    (*ph_orig)[idx].nearest_block_index=-1;
                    (*ph_orig)[idx].weight=0; //now making the absorbed UNABSORBED_CS_PHOTON photons become null photons

                }
                else
                {
                    //all rebinned photns have been saved so just treat the rest of the phootn array as null photons
                    (*ph_orig)[idx].weight=0;
                    (*ph_orig)[idx].nearest_block_index=-1;
                    //j++;
                }
            }
            
            
        }
        else if ((*ph_orig)[idx].type != CS_POOL_PHOTON)
        {
            //this is a realloc photon that has to be set to null
            (*ph_orig)[idx].type = COMPTONIZED_PHOTON;
            (*ph_orig)[idx].weight=0;
            (*ph_orig)[idx].nearest_block_index=-1;

        }
        
        //if ((*ph_orig)[idx].weight==0)
        //{
            //if the bin had no photons in it, the rebinned photon is effectively null
        //    j++;
        //}
        
    }
    
    //make sure that all the rebinned photons have been saved
    if (count<num_bins)
    {
        fprintf(fPtr, "There was an issue where MCRaT was not able to save all of the rebinned photons\n");
        printf("TThere was an issue where MCRaT was not able to save all of the rebinned photons\n");
        fflush(fPtr);
        exit(1);
    }
    //fprintf(fPtr, "i count after first loop %d\n", idx+1);
    //make sure we look at whole array of photons to see hwo many null photons we have
    /*
    for (i=idx+1;i < *num_ph; i++)
    {
        if ((*ph_orig)[i].weight==0)
        {
            //if the bin had no photons in it, the rebinned photon is effectively null
            j++;
            fprintf(fPtr, "i count in the if %d\n", i);
        }
        
        //i++;
    }
     */
    //fprintf(fPtr, "i count after second loop %d\n", i);
    
    
    
    int null_ph_count=0;
    int null_ph_count_1=0;
    //int null_ph_count_2=0;
#pragma omp parallel for num_threads(num_thread) reduction(+:null_ph_count)
    for (i=0;i<*num_ph;i++)
    {
        if ((*ph_orig)[i].weight == 0)
        {
            null_ph_count++;
            //fprintf(fPtr, "%d \n", null_ph_count);
            
        }
        /*
        if ((*ph_orig)[i].type == INJECTED_PHOTON)
        {
            null_ph_count_1++;
        }
          */
        //fprintf(fPtr, "%d %c %e %e %e\n", i, (*ph_orig)[i].type, (*ph_orig)[i].weight, (*ph_orig)[i].p0, (*ph_orig)[i].s0 );
        
    }
    
    //*num_cyclosynch_ph_emit=0;
    
    //int temporary;
    *scatt_cyclosynch_num_ph=num_bins-num_null_rebin_ph;
    //temporary=num_bins+synch_photon_count-num_null_rebin_ph;
    *num_cyclosynch_ph_emit=num_bins+synch_photon_count-num_null_rebin_ph; //include the emitted synch photons and exclude any of those that are null
    *num_null_ph=null_ph_count; //was using j before but i have no idea why its not counting correctly
        
    //fprintf(fPtr, "orig null_ph: %d Calc num_ph: %d counted null_ph: %d forloop null_ph: %d, num_inj: %d num_null_rebin_ph: %d old scatt_cyclosynch_num_ph: %d new scatt_cyclosynch_num_ph: %d\n", *num_null_ph, (*num_ph), j, null_ph_count, null_ph_count_1, num_null_rebin_ph, *scatt_cyclosynch_num_ph, num_bins-num_null_rebin_ph  );
    
    //fprintf(fPtr, "at end of rebin f(x)  %d,  %d, %d\n", null_ph_count, num_bins+synch_photon_count-num_null_rebin_ph, num_bins-num_null_rebin_ph);
    //fflush(fPtr);
    /*
    count=0;
    for (i=0;i<*num_ph;i++)
    {
        //fprintf(fPtr, "%d %c %e %e\n", i, (*ph_orig)[i].type, (*ph_orig)[i].weight, (*ph_orig)[i].p0 );
        //fflush(fPtr);
        
        if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
        {
            count++;
        }
    }
    fprintf(fPtr, "at end of rebin f(x): count is: %d and scatt_cyclosynch_num_ph is %d\n", count, *scatt_cyclosynch_num_ph );
    */
    
    ////gsl_histogram_fprintf (stdout, h, "%g", "%g");
    gsl_histogram_free (h); gsl_histogram_free(h_theta_space); gsl_histogram2d_free (h_energy_theta);
    gsl_histogram2d_pdf_free (pdf_phi_theta);
    gsl_histogram2d_free (h_phi_theta);
    free(rebin_ph);
    //free( synch_comp_photon_idx);
 
    return num_null_rebin_ph; //num_bins-num_null_rebin_ph;
}

int rebin2dCyclosynchCompPhotons(struct photon **ph_orig, int *num_ph,  int *num_null_ph, int *num_cyclosynch_ph_emit, int *scatt_cyclosynch_num_ph, double **all_time_steps, int **sorted_indexes, int max_photons, double thread_theta_min, double thread_theta_max , gsl_rng * rand, FILE *fPtr)
{
    int i=0, j=0, k=0, count=0, count_x=0, count_y=0, count_c_ph=0, end_count=(*scatt_cyclosynch_num_ph), idx=0, num_thread=1;
    #if defined(_OPENMP)
    num_thread=omp_get_num_threads();
    #endif
    int synch_comp_photon_count=0, synch_photon_count=0, num_avg=12, num_bins=(CYCLOSYNCHROTRON_REBIN_E_PERC)*max_photons; //some factor of the max number of photons that is specified in the mc.par file, num bins is also in test function
    double dtheta_bin=CYCLOSYNCHROTRON_REBIN_ANG*M_PI/180; //the size of the bin that we want to produce for spatial binning in theta
    int num_bins_theta=(thread_theta_max-thread_theta_min)/dtheta_bin;//try this many bins such that we have 0.5 degree resolution, can also try to do adaptive binning with constant SNR
    double avg_values[12]={0}; //number of averages that'll be taken is given by num_avg in above line
    double p0_min=DBL_MAX, p0_max=0, log_p0_min=0, log_p0_max=0;//look at p0 of photons not by frequency since its just nu=p0*C_LIGHT/PL_CONST
    double rand1=0, rand2=0, phi=0, theta=0;
    double min_range=0, max_range=0, min_range_theta=0, max_range_theta=0, energy=0;
    double ph_r=0, ph_theta=0, temp_theta_max=0, temp_theta_min=DBL_MAX;
    //int *synch_comp_photon_idx=NULL; make this an array b/c had issue with deallocating this memory for some reason
    int synch_comp_photon_idx[*scatt_cyclosynch_num_ph];
    //struct photon *rebin_ph=malloc(num_bins* sizeof (struct photon ));
    int num_null_rebin_ph=0, num_in_bin=0;
    struct photon *tmp=NULL;
    double *tmp_double=NULL;
    int *tmp_int=NULL;
    double count_weight=0;

    fprintf(fPtr, "In the rebin func; num_threads %d scatt_cyclosynch_num_ph %d, num_ph %d\n", num_thread, (*scatt_cyclosynch_num_ph), *num_ph);
    fflush(fPtr);
    
    int min_idx=0, max_idx=0;
    count=0;
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
        {
            //see if the photon's nu is larger than nu_max or smaller than nu_min
            if (((*ph_orig)[i].p0< p0_min))
            {
                //dont include any absorbed UNABSORBED_CS_PHOTON photons that have negative P0 values
                p0_min= (*ph_orig)[i].p0;
                min_idx=i;
                //fprintf(fPtr, "new p0 min %e\n", (p0_min) );
            }
            
            if ((*ph_orig)[i].p0> p0_max)
            {
                p0_max= (*ph_orig)[i].p0;
                max_idx=i;
                //fprintf(fPtr, "new p0 max %e\n", (p0_max) );
            }
            
            //look at min and max theta of photons
            ph_r=pow(((*ph_orig)[i].r0)*((*ph_orig)[i].r0) + ((*ph_orig)[i].r1)*((*ph_orig)[i].r1) + ((*ph_orig)[i].r2)*((*ph_orig)[i].r2),0.5);
            ph_theta=acos(((*ph_orig)[i].r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            
            if (ph_theta > temp_theta_max )
            {
                temp_theta_max=ph_theta;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_theta<temp_theta_min)
            {
                temp_theta_min=ph_theta;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }

            
            // also save the index of these photons because they wil become null later on
            synch_comp_photon_idx[count]=i;
            //fprintf(fPtr, "Save index %d\n", i );
            count++;
            
            if ((*ph_orig)[i].type == COMPTONIZED_PHOTON)
            {
                //keep track of the number of COMPTONIZED_PHOTON photons so we can know if the array needs to be increased in size, also take num_null_ph into account in doing this
                count_c_ph+=1;
            }
        }
        else if (((*ph_orig)[i].type == CS_POOL_PHOTON) && ((*ph_orig)[i].weight != 0))
        {
            synch_photon_count++;
        }
    }
    
    
    num_bins_theta=1+(temp_theta_max-temp_theta_min)/dtheta_bin;
    
    fprintf(fPtr, "Rebin: min, max (keV): %e %e log p0 min, max: %e %e idx: %d %d\n", p0_min*C_LIGHT/1.6e-9,p0_max*C_LIGHT/1.6e-9 , log10(p0_min), log10(p0_max), min_idx, max_idx );
    fprintf(fPtr, "Rebin: min, max (theta in deg): %e %e number of bins %d count: %d\n", temp_theta_min*180/M_PI, temp_theta_max*180/M_PI, num_bins_theta, count );
    fflush(fPtr);
    
    if (count != end_count)
    {
        end_count=count; //need this for some reason idk why end_count gets off by 1 compared to what it should be
        fprintf(fPtr, "Rebin: not equal to end_count therefore resetting count to be: %d\n", count );
        fflush(fPtr);
    }
    
    if (num_bins_theta*num_bins>=max_photons)
    {
        fprintf(fPtr, "The number of rebinned photons, %d, is larger than max_photons %d and will not rebin efficiently. Adjust the parameters such that the number of bins in theta and energy are less than the number of photons that will lead to rebinning.\n",  num_bins_theta*num_bins, max_photons);
        fflush(fPtr);
        
        printf("Rebin: min, max (theta in deg): %e %e number of bins %d count: %d\n", temp_theta_min*180/M_PI, temp_theta_max*180/M_PI, num_bins_theta, count );
        printf( "In angle range: %e-%e: The number of rebinned photons, %d, is larger than max_photons %d and will not rebin efficiently. Adjust the parameters such that the number of bins in theta and energy are less than the number of photons that will lead to rebinning.\n",  thread_theta_min*180/M_PI, thread_theta_max*180/M_PI, num_bins_theta*num_bins, max_photons);
        exit(1);

    }

        
    struct photon *rebin_ph=malloc(num_bins*num_bins_theta* sizeof (struct photon ));
    struct photon *synch_ph=malloc(synch_photon_count* sizeof (struct photon ));
    int synch_photon_idx[synch_photon_count];
    
    gsl_histogram2d * h_energy_theta = gsl_histogram2d_alloc (num_bins, num_bins_theta); //x is for energy  and y is for spatial theta, goes from 0 to pi
    gsl_histogram2d_set_ranges_uniform (h_energy_theta, log10(p0_min), log10(p0_max*(1+1e-6)), temp_theta_min, temp_theta_max+dtheta_bin);

    //populate histogram for photons with nu that falss within the proper histogram bin
    //may not need this loop, can just check if the photon nu falls within the bin edges and do averages etc within next loop
    count=0;
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
        {
            
            ph_r=pow(((*ph_orig)[i].r0)*((*ph_orig)[i].r0) + ((*ph_orig)[i].r1)*((*ph_orig)[i].r1) + ((*ph_orig)[i].r2)*((*ph_orig)[i].r2),0.5);
            ph_theta=acos(((*ph_orig)[i].r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            
            gsl_histogram2d_increment(h_energy_theta, log10((*ph_orig)[i].p0), ph_theta);
            
            count_weight+=(*ph_orig)[i].weight;

        }
        
        if (((*ph_orig)[i].type == CS_POOL_PHOTON) && ((*ph_orig)[i].weight != 0))
        {
            //save the sych photons here because they may get written over later and corrupted
            (synch_ph+count)->p0=(*ph_orig)[i].p0;
            (synch_ph+count)->p1=(*ph_orig)[i].p1;
            (synch_ph+count)->p2=(*ph_orig)[i].p2;
            (synch_ph+count)->p3=(*ph_orig)[i].p3;
            (synch_ph+count)->comv_p0=(*ph_orig)[i].comv_p0;
            (synch_ph+count)->comv_p1=(*ph_orig)[i].comv_p1;
            (synch_ph+count)->comv_p2=(*ph_orig)[i].comv_p2;
            (synch_ph+count)->comv_p3=(*ph_orig)[i].comv_p3;
            (synch_ph+count)->r0=(*ph_orig)[i].r0;
            (synch_ph+count)->r1= (*ph_orig)[i].r1;
            (synch_ph+count)->r2=(*ph_orig)[i].r2; //y coordinate in flash becomes z coordinate in MCRaT
            (synch_ph+count)->s0=(*ph_orig)[i].s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
            (synch_ph+count)->s1=(*ph_orig)[i].s1;
            (synch_ph+count)->s2=(*ph_orig)[i].s2;
            (synch_ph+count)->s3=(*ph_orig)[i].s3;
            (synch_ph+count)->num_scatt=(*ph_orig)[i].num_scatt;
            (synch_ph+count)->weight=(*ph_orig)[i].weight;
            (synch_ph+count)->nearest_block_index=(*ph_orig)[i].nearest_block_index; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
            synch_photon_idx[count]=i;
            count++;
        }
         
    }
    
    //fprintf(fPtr, "counted_weight 1 %e\n", count_weight);
    //fflush(fPtr);
    //count_weight=0;
    //gsl_histogram2d_fprintf(fPtr, h_energy_theta, "%g", "%g");
        
    for (count_x=0;count_x<num_bins;count_x++)
    {
        for (count_y=0;count_y<num_bins_theta;count_y++)
        {
            count=count_x*num_bins_theta+count_y;
            gsl_histogram2d_get_xrange(h_energy_theta, count_x, &min_range, &max_range);
            gsl_histogram2d_get_yrange(h_energy_theta, count_y, &min_range_theta, &max_range_theta);
            
            num_in_bin=gsl_histogram2d_get(h_energy_theta, count_x, count_y);
            
            if (num_in_bin > 1)
            {
                
                rand1=gsl_rng_uniform(rand)*num_in_bin;//random photon to choose theta and phi of 4 mometum for rebinned photon
                
                for (j=0;j<num_avg;j++)
                {
                    avg_values[j]=0;
                }
                
                //loop over the number of photons
                j=0; //to compare to rand1 for choosing the photon
                for (i=0;i<*num_ph;i++)
                {
                    if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON)) && ((*ph_orig)[i].p0 > 0))
                    {
                        ph_r=pow(((*ph_orig)[i].r0)*((*ph_orig)[i].r0) + ((*ph_orig)[i].r1)*((*ph_orig)[i].r1) + ((*ph_orig)[i].r2)*((*ph_orig)[i].r2),0.5);
                        ph_theta=acos(((*ph_orig)[i].r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians

                        //if the photon nu falls in the count bin of the nu histogram then add it to the phi_theta 2d hist
                        if ((log10((*ph_orig)[i].p0)< max_range  ) && (log10((*ph_orig)[i].p0)>=min_range) && (ph_theta < max_range_theta) && (ph_theta >= min_range_theta))
                        {
                            //not doing average values, choosing random photon instead
                            avg_values[0] += ph_r*(*ph_orig)[i].weight; // doing r, theta averages in space
                            avg_values[1] += ph_theta*(*ph_orig)[i].weight;
                            avg_values[2] += ((atan((*ph_orig)[i].p2/((*ph_orig)[i].p1))*180/M_PI)-(atan(((*ph_orig)[i].r1)/ ((*ph_orig)[i].r0))*180/M_PI))*(*ph_orig)[i].weight;// look at delta \phi between the 4 mometum and its location
                            avg_values[3] += (*ph_orig)[i].s0*(*ph_orig)[i].weight;
                            avg_values[4] += (*ph_orig)[i].s1*(*ph_orig)[i].weight;
                            avg_values[5] += (*ph_orig)[i].s2*(*ph_orig)[i].weight;
                            avg_values[6] += (*ph_orig)[i].s3*(*ph_orig)[i].weight;
                            avg_values[7] += (*ph_orig)[i].num_scatt*(*ph_orig)[i].weight;
                            avg_values[8] += (*ph_orig)[i].weight;
                            
                            //get theta and phi of random photon
                            {
                                avg_values[9] += fmod(atan2((*ph_orig)[i].p2,((*ph_orig)[i].p1))*180/M_PI + 360.0,360.0) *(*ph_orig)[i].weight;
                                avg_values[10] += (180/M_PI)*acos(((*ph_orig)[i].p3)/((*ph_orig)[i].p0))*(*ph_orig)[i].weight;

                            }
                            
                            avg_values[11] +=(*ph_orig)[i].p0*(*ph_orig)[i].weight;
                            
                            j++;
                            
                        }
                    }
                }
                                               
                
                energy=avg_values[11]/avg_values[8];
                
                
                phi=avg_values[9]/avg_values[8];
                theta=avg_values[10]/avg_values[8];
                
                (rebin_ph+count)->type = COMPTONIZED_PHOTON;
                                
                
                (rebin_ph+count)->p0=energy;
                (rebin_ph+count)->p1=energy*sin(theta*M_PI/180)*cos(phi*M_PI/180);
                (rebin_ph+count)->p2=energy*sin(theta*M_PI/180)*sin(phi*M_PI/180);
                (rebin_ph+count)->p3=energy*cos(theta*M_PI/180);
                (rebin_ph+count)->comv_p0=0;
                (rebin_ph+count)->comv_p1=0;
                (rebin_ph+count)->comv_p2=0;
                (rebin_ph+count)->comv_p3=0;
                
                //calculate the rebinned photon's positional phi as a displacement from its 4-momentum phi direction
                rand1=(M_PI/180)*(phi-avg_values[2]/avg_values[8]);
                
                
                (rebin_ph+count)->r0= (avg_values[0]/avg_values[8])*sin(avg_values[1]/avg_values[8])*cos(rand1); //avg_values[0]/avg_values[8]; now do avg r * avg theta * random phi
                (rebin_ph+count)->r1= (avg_values[0]/avg_values[8])*sin(avg_values[1]/avg_values[8])*sin(rand1); //avg_values[1]/avg_values[8];
                (rebin_ph+count)->r2= (avg_values[0]/avg_values[8])*cos(avg_values[1]/avg_values[8]); //avg_values[2]/avg_values[8];

                (rebin_ph+count)->s0=avg_values[3]/avg_values[8]; // stokes parameterized are normalized such that I always =1
                (rebin_ph+count)->s1=avg_values[4]/avg_values[8];
                (rebin_ph+count)->s2=avg_values[5]/avg_values[8];
                (rebin_ph+count)->s3=avg_values[6]/avg_values[8];
                (rebin_ph+count)->num_scatt=avg_values[7]/avg_values[8];
                (rebin_ph+count)->weight=avg_values[8];
                (rebin_ph+count)->nearest_block_index=0; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
                
                //fprintf(fPtr, "bin %e-%e, %e-%e has %e photons: Theta of averages photon is: %e\n",pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, min_range_theta, max_range_theta, gsl_histogram2d_get(h_energy_theta, count_x, count_y), ph_theta);
                //fflush(fPtr);
                
                count_weight+=(rebin_ph+count)->weight;
                                
            }
            else if (num_in_bin == 1)
            {
                //fprintf(fPtr, "bin %e-%e has %e photons\n", pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, 1.0);

                //for thr case of just 1 hoton being in the bin just set the rebinned photon to the one photons parameters
                for (i=0;i<*num_ph;i++)
                {
                    if (((*ph_orig)[i].weight != 0) && (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON))&& ((*ph_orig)[i].p0 > 0))
                    {
                        ph_r=pow(((*ph_orig)[i].r0)*((*ph_orig)[i].r0) + ((*ph_orig)[i].r1)*((*ph_orig)[i].r1) + ((*ph_orig)[i].r2)*((*ph_orig)[i].r2),0.5);
                        ph_theta=acos(((*ph_orig)[i].r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians

                        //if the photon nu falls in the count bin of the nu histogram then add it to the phi_theta 2d hist
                        if ((log10((*ph_orig)[i].p0)< max_range  ) && (log10((*ph_orig)[i].p0)>=min_range) && (ph_theta < max_range_theta) && (ph_theta >= min_range_theta))
                        {
                            (rebin_ph+count)->p0=(*ph_orig)[i].p0;
                            (rebin_ph+count)->p1=(*ph_orig)[i].p1;
                            (rebin_ph+count)->p2=(*ph_orig)[i].p2;
                            (rebin_ph+count)->p3=(*ph_orig)[i].p3;
                            (rebin_ph+count)->comv_p0=(*ph_orig)[i].comv_p0;
                            (rebin_ph+count)->comv_p1=(*ph_orig)[i].comv_p1;
                            (rebin_ph+count)->comv_p2=(*ph_orig)[i].comv_p2;
                            (rebin_ph+count)->comv_p3=(*ph_orig)[i].comv_p3;
                            (rebin_ph+count)->r0=(*ph_orig)[i].r0;
                            (rebin_ph+count)->r1= (*ph_orig)[i].r1;
                            (rebin_ph+count)->r2=(*ph_orig)[i].r2; //y coordinate in flash becomes z coordinate in MCRaT
                            (rebin_ph+count)->s0=(*ph_orig)[i].s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                            (rebin_ph+count)->s1=(*ph_orig)[i].s1;
                            (rebin_ph+count)->s2=(*ph_orig)[i].s2;
                            (rebin_ph+count)->s3=(*ph_orig)[i].s3;
                            (rebin_ph+count)->num_scatt=(*ph_orig)[i].num_scatt;
                            (rebin_ph+count)->weight=(*ph_orig)[i].weight;
                            (rebin_ph+count)->nearest_block_index=(*ph_orig)[i].nearest_block_index; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
                            
                            (rebin_ph+count)->type = COMPTONIZED_PHOTON;
                            
                            count_weight+=(rebin_ph+count)->weight;

                            i=*num_ph;
                        }
                    }
                }
            }
            else
            {
                //fprintf(fPtr, "Rebinned Photon is a null photon because there are no photons in this energy bin.\n");
                (rebin_ph+count)->type = COMPTONIZED_PHOTON;
                
                energy=pow(10,0.5*(max_range+min_range));
                
                //fprintf(fPtr, "bin %e-%e has %e photons\n", pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, 0.0);

                
                (rebin_ph+count)->p0=energy;
                (rebin_ph+count)->p1=0;
                (rebin_ph+count)->p2=0;
                (rebin_ph+count)->p3=0;
                (rebin_ph+count)->comv_p0=0;
                (rebin_ph+count)->comv_p1=0;
                (rebin_ph+count)->comv_p2=0;
                (rebin_ph+count)->comv_p3=0;
                (rebin_ph+count)->r0=0;
                (rebin_ph+count)->r1= 0;
                (rebin_ph+count)->r2=0;
                (rebin_ph+count)->s0=1; // stokes parameterized are normalized such that I always =1
                (rebin_ph+count)->s1=0;
                (rebin_ph+count)->s2=0;
                (rebin_ph+count)->s3=0;
                (rebin_ph+count)->num_scatt=0;
                (rebin_ph+count)->weight=0;
                (rebin_ph+count)->nearest_block_index=-1; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
                
                count_weight+=(rebin_ph+count)->weight;


            }
            
            
            
        }
    }
    
    //fprintf(fPtr, "counted_weight 2 %e\n", count_weight);
    //fflush(fPtr);
    
    //exit(0);
    
    
    if ((count_c_ph+(*num_null_ph))<num_bins*num_bins_theta)
    {
        //need to expand the array
        //if the totoal number of photons to be emitted is larger than the number of null phtons curently in the array, then have to grow the array
        //need to realloc memory to hold the old photon info and the new emitted photon's info
        //before was doing ((*num_ph)+num_bins-count_c_ph-(*num_null_ph) ) but instead did the beow since that creatd extr space for null photons to be filled in later on
        
        /*
         fprintf(fPtr, "\nBefore Rebin: \n");
         fflush(fPtr);

        for (i=0;i<*num_ph;i++)
        {
            fprintf(fPtr, "%d %c %e %e %e %d\n", i, (*ph_orig)[i].type, (*ph_orig)[i].weight, (*ph_orig)[i].p0, (*ph_orig)[i].s0, (*ph_orig)[i].nearest_block_index );

        }
         fprintf(fPtr, "Before Rebin: \n\n");
         fflush(fPtr);
         */

        
        //fprintf(fPtr, "Rebin: Allocating %d space\n", ((*num_ph)+num_bins*num_bins_theta-count_c_ph+(*num_null_ph) )); //befoe was dong
        //fflush(fPtr);
        tmp=realloc(*ph_orig, ((*num_ph)+num_bins*num_bins_theta-count_c_ph+(*num_null_ph))* sizeof (struct photon )); //may have to look into directly doubling (or *1.5) number of photons each time we need to allocate more memory, can do after looking at profiling for "just enough" memory method
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
            exit(1);
        }
        //fprintf(fPtr, "Rebin: Allocating %d space good in 1st realloc\n", ((*num_ph)+num_bins-count_c_ph+(*num_null_ph) ));
        //fflush(fPtr);

        
        //also expand memory of other arrays
        tmp_double=realloc(*all_time_steps, ((*num_ph)+num_bins*num_bins_theta-count_c_ph+(*num_null_ph))*sizeof(double));
        if (tmp_double!=NULL)
        {
            *all_time_steps=tmp_double;
        }
        else
        {
            printf("Error with reallocating space to hold data about each photon's time step until an interaction occurs\n");
            exit(1);

        }
        
        //fprintf(fPtr, "Rebin: Allocating %d space good in 2nd realloc\n", ((*num_ph)+num_bins-count_c_ph+(*num_null_ph) ));
        //fflush(fPtr);

        
        tmp_int=realloc(*sorted_indexes, ((*num_ph)+num_bins*num_bins_theta-count_c_ph+(*num_null_ph))*sizeof(int));
        if (tmp_int!=NULL)
        {
            *sorted_indexes=tmp_int;
        }
        else
        {
            printf("Error with reallocating space to hold data about the order in which each photon would have an interaction\n");
            exit(1);

        }

        
        *num_ph=( *num_ph)+(num_bins*num_bins_theta-count_c_ph+(*num_null_ph));
        end_count=(*scatt_cyclosynch_num_ph)+num_bins*num_bins_theta-count_c_ph+(*num_null_ph);
        
        //go through all the phtons and ID where the synchrotron ones are and reset them to what they were originally
        count=0;
        for (i=0;i<synch_photon_count;i++)
        {
            //if ((*ph_orig)[i].type == CS_POOL_PHOTON)
            {
                idx=synch_photon_idx[i];
                (*ph_orig)[idx].p0=(synch_ph+count)->p0;
                (*ph_orig)[idx].p1=(synch_ph+count)->p1;
                (*ph_orig)[idx].p2=(synch_ph+count)->p2;
                (*ph_orig)[idx].p3=(synch_ph+count)->p3;
                (*ph_orig)[idx].comv_p0=(synch_ph+count)->comv_p0;
                (*ph_orig)[idx].comv_p1=(synch_ph+count)->comv_p1;
                (*ph_orig)[idx].comv_p2=(synch_ph+count)->comv_p2;
                (*ph_orig)[idx].comv_p3=(synch_ph+count)->comv_p3;
                (*ph_orig)[idx].r0=(synch_ph+count)->r0;
                (*ph_orig)[idx].r1=(synch_ph+count)->r1;
                (*ph_orig)[idx].r2=(synch_ph+count)->r2; //y coordinate in flash becomes z coordinate in MCRaT
                (*ph_orig)[idx].s0=(synch_ph+count)->s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                (*ph_orig)[idx].s1=(synch_ph+count)->s1;
                (*ph_orig)[idx].s2=(synch_ph+count)->s2;
                (*ph_orig)[idx].s3=(synch_ph+count)->s3;
                (*ph_orig)[idx].num_scatt=(synch_ph+count)->num_scatt;
                (*ph_orig)[idx].weight=(synch_ph+count)->weight;
                (*ph_orig)[idx].nearest_block_index=(synch_ph+count)->nearest_block_index;
                //(*ph_orig)[idx].type = CS_POOL_PHOTON;
                
                count++;
            }
        }
        
        //go through the newly added phootns in the array and make sure that there are no CS_POOL_PHOTON type photons there that can cause issues later on
        for (i=( *num_ph)-(num_bins*num_bins_theta-count_c_ph+(*num_null_ph));i<(*num_ph);i++)
        {
            if ((*ph_orig)[i].type == CS_POOL_PHOTON)
            {
                (*ph_orig)[i].type = COMPTONIZED_PHOTON;
            }
        }
        
        /*
         fprintf(fPtr, "\nAfter Rebin: \n");
         fflush(fPtr);

        for (i=0;i<*num_ph;i++)
        {
            fprintf(fPtr, "%d %c %e %e %e %d\n", i, (*ph_orig)[i].type, (*ph_orig)[i].weight, (*ph_orig)[i].p0, (*ph_orig)[i].s0, (*ph_orig)[i].nearest_block_index );

        }
         fprintf(fPtr, "After Rebin: \n\n");
         fflush(fPtr);
         */

        
        
    }
    
    //go through and assign the rebinned photons to the COMPTONIZED_PHOTON phtoons and make all UNABSORBED_CS_PHOTON photons become "absorbed"
    j=0;
    count=0;
    i=0;
    for (i=0;i<end_count;i++)
    {
        if (i<(*scatt_cyclosynch_num_ph))
        {
            //the photon idx can be found from the original array
            idx=synch_comp_photon_idx[i];
        }
        else
        {
            //enter this if we realloc the arrays
            idx=i-(*scatt_cyclosynch_num_ph)+( *num_ph)-(num_bins*num_bins_theta-count_c_ph+(*num_null_ph));//go to the end of the old num_ph value and start setting things to null photons
        }

        if (((*ph_orig)[idx].type == UNABSORBED_CS_PHOTON) || ((*ph_orig)[idx].type == COMPTONIZED_PHOTON))
        {
            if (count<num_bins*num_bins_theta)
            {
                //keep saving the rebinned photns
                (*ph_orig)[idx].p0=(rebin_ph+count)->p0;
                (*ph_orig)[idx].p1=(rebin_ph+count)->p1;
                (*ph_orig)[idx].p2=(rebin_ph+count)->p2;
                (*ph_orig)[idx].p3=(rebin_ph+count)->p3;
                (*ph_orig)[idx].comv_p0=(rebin_ph+count)->comv_p0;
                (*ph_orig)[idx].comv_p1=(rebin_ph+count)->comv_p1;
                (*ph_orig)[idx].comv_p2=(rebin_ph+count)->comv_p2;
                (*ph_orig)[idx].comv_p3=(rebin_ph+count)->comv_p3;
                (*ph_orig)[idx].r0=(rebin_ph+count)->r0;
                (*ph_orig)[idx].r1=(rebin_ph+count)->r1;
                (*ph_orig)[idx].r2=(rebin_ph+count)->r2; //y coordinate in flash becomes z coordinate in MCRaT
                (*ph_orig)[idx].s0=(rebin_ph+count)->s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                (*ph_orig)[idx].s1=(rebin_ph+count)->s1;
                (*ph_orig)[idx].s2=(rebin_ph+count)->s2;
                (*ph_orig)[idx].s3=(rebin_ph+count)->s3;
                (*ph_orig)[idx].num_scatt=(rebin_ph+count)->num_scatt;
                (*ph_orig)[idx].weight=(rebin_ph+count)->weight;
                (*ph_orig)[idx].nearest_block_index=(rebin_ph+count)->nearest_block_index;
                (*ph_orig)[idx].type = COMPTONIZED_PHOTON;
                                
                if ((rebin_ph+count)->weight==0)
                {
                    //if the bin had no photons in it, the rebinned photon is effectively null
                    //j++;
                    num_null_rebin_ph++;
                }
                
                count++;
                
            }
            else
            {
                if ((*ph_orig)[idx].type == UNABSORBED_CS_PHOTON)
                {
                    //if the photon is a compton scatered synch photon from a past frame treat it as though it has been absorbed
                    (*ph_orig)[idx].p0=-1; //set its energy negative so we know for later analysis that it can't be used and its been "absorbed", this makes it still get saves in the hdf5 files
                    (*ph_orig)[idx].nearest_block_index=-1;
                    (*ph_orig)[idx].weight=0; //now making the absorbed UNABSORBED_CS_PHOTON photons become null photons
                    //num_null_rebin_ph++;
                }
                else
                {
                    //all rebinned photns have been saved so just treat the rest of the phootn array as null photons
                    (*ph_orig)[idx].weight=0;
                    (*ph_orig)[idx].nearest_block_index=-1;
                    //j++;
                    //num_null_rebin_ph++;
                }
            }
            
            
        }
        else if ((*ph_orig)[idx].type != CS_POOL_PHOTON)
        {
            //this is a realloc photon that has to be set to null
            (*ph_orig)[idx].type = COMPTONIZED_PHOTON;
            (*ph_orig)[idx].weight=0;
            (*ph_orig)[idx].nearest_block_index=-1;
            //num_null_rebin_ph++;
        }
        
        if ((*ph_orig)[idx].weight==0)
        {
            //if the bin had no photons in it, the rebinned photon is effectively null
            j++;
        }
        
    }
    
    //make sure that all the rebinned photons have been saved
    if (count<num_bins*num_bins_theta)
    {
        fprintf(fPtr, "There was an issue where MCRaT was not able to save all of the rebinned photons\n");
        printf("TThere was an issue where MCRaT was not able to save all of the rebinned photons\n");
        fflush(fPtr);
        exit(1);
    }
    
    /*
     fprintf(fPtr, "\nPost Rebin Fill-in: \n");
     fflush(fPtr);

    for (i=0;i<*num_ph;i++)
    {
        fprintf(fPtr, "%d %c %e %e %e %d\n", i, (*ph_orig)[i].type, (*ph_orig)[i].weight, (*ph_orig)[i].p0, (*ph_orig)[i].s0, (*ph_orig)[i].nearest_block_index );

    }
     fprintf(fPtr, "Post Rebin Fill-in: \n\n");
     fflush(fPtr);

     */

    
    
    *scatt_cyclosynch_num_ph=num_bins*num_bins_theta-num_null_rebin_ph;
    *num_cyclosynch_ph_emit=num_bins*num_bins_theta+synch_photon_count-num_null_rebin_ph; //include the emitted synch photons and exclude any of those that are null
    *num_null_ph=j; //was using j before but i have no idea why its not counting correctly
        
    //fprintf(fPtr, "orig null_ph: %d Calc num_ph: %d counted null_ph: %d  num_null_rebin_ph: %d old scatt_cyclosynch_num_ph: %d new scatt_cyclosynch_num_ph: %d\n", *num_null_ph, (*num_ph), j, num_null_rebin_ph, *scatt_cyclosynch_num_ph, num_bins*num_bins_theta-num_null_rebin_ph  );
    //fflush(fPtr);
    
    //exit(2);

    
    gsl_histogram2d_free (h_energy_theta);
    free(rebin_ph);
    free(synch_ph);
    
    return num_null_rebin_ph;
}


int photonEmitCyclosynch(struct photon **ph_orig, int *num_ph, int *num_null_ph, double **all_time_steps, int **sorted_indexes, double r_inj, double ph_weight, int maximum_photons, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng *rand, int inject_single_switch, int scatt_ph_index, FILE *fPtr)
{
    double rmin=0, rmax=0, max_photons=CYCLOSYNCHROTRON_REBIN_E_PERC*maximum_photons; //have 10% as default, can change later need to figure out how many photons across simulations I want emitted
    double ph_weight_adjusted=0, position_phi=0;
    double dimlesstheta=0, nu_c=0, el_dens=0, error=0, ph_dens_calc=0, max_jnu=0;
    double r_grid_innercorner=0, r_grid_outercorner=0, theta_grid_innercorner=0, theta_grid_outercorner=0;
    double el_p[4], ph_p_comv[4];
    double params[3];
    double fr_dum=0.0, y_dum=0.0, yfr_dum=0.0, com_v_phi=0, com_v_theta=0, position_rand=0, position2_rand=0, position3_rand=0, cartesian_position_rand_array[3];
    double *p_comv=NULL, *boost=NULL, *l_boost=NULL; //pointers to hold comov 4 monetum, the fluid vlocity, and the photon 4 momentum in the lab frame
    int status;
    int block_cnt=0, i=0, j=0, k=0, null_ph_count=0, *ph_dens=NULL, ph_tot=0, net_ph=0, min_photons=1;
    int *null_ph_indexes=NULL;
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    int count_null_indexes=0, idx=0;
    struct photon *ph_emit=NULL; //pointer to array of structs that will hold emitted photon info
    struct photon *tmp=NULL;
    double *tmp_double=NULL;
    int *tmp_int=NULL, n_pool=0;
    
    //fprintf(fPtr, "IN EMIT SYNCH FUNCTION; num_threads %d\n", num_thread);
    //fprintf(fPtr, "BEFORE Original number of photons: %d Null photons %d\n", (*num_ph), null_ph_count, ph_tot);
    //fflush(fPtr);
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (10000);
    
    gsl_function F;
    F.function = &blackbody_ph_spect; //&jnu_ph_spect;
    
    if (inject_single_switch == 0)
    {
        rmin=calcCyclosynchRLimits(hydro_data->scatt_frame_number,  hydro_data->inj_frame_number, hydro_data->fps,  r_inj, "min");
        rmax=calcCyclosynchRLimits(hydro_data->scatt_frame_number,  hydro_data->inj_frame_number, hydro_data->fps,  r_inj, "max");
        
        fprintf(fPtr, "rmin %e rmax %e, theta min/max: %e %e\n", rmin, rmax, theta_min, theta_max);
        #pragma omp parallel for num_threads(num_thread) reduction(+:block_cnt)
        for(i=0;i<hydro_data->num_elements;i++)
        {
            //look at all boxes in width delta r=c/fps and within angles we are interested in NEED TO IMPLEMENT
            #if DIMENSIONS == 3
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]-0.5*(hydro_data->r2_size)[i]);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]+0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
            #endif

            if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  < rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner < theta_max))
            {
                block_cnt+=1;
            }
            
        }
        
        fprintf(fPtr, "MCRaT has chosen %d hydro elements that it will emit cyclosynchrotron photons into.\n", block_cnt);
        fflush(fPtr);
        
        //min_photons=block_cnt;//do this so we have at least one synch photon in each block that meets the radius requiements, need to double check to see if this changes things or not (I dont expect it to), testing shows that this may not work as block_cnt can be larger than max_photons
        
        if (block_cnt==0)
        {
            min_photons=block_cnt; //this is for the case of there being no blocks near photons, probably b/c photons are off the hydro grid, and lets the program progress past the while loop below
        }
        
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
            
            for (i=0;i< hydro_data->num_elements;i++)
            {
                //printf("%d\n",i);
                //printf("%e, %e, %e, %e, %e, %e\n", *(r+i),(r_inj - C_LIGHT/fps), (r_inj + C_LIGHT/fps), *(theta+i) , theta_max, theta_min);
                #if DIMENSIONS == 3
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]-0.5*(hydro_data->r2_size)[i]);
                        hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]+0.5*(hydro_data->r2_size)[i]);
                #else
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
                #endif

                if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  < rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner < theta_max))
                {
                    //set parameters fro integration fo phtoons spectrum
                    el_dens= ((hydro_data->dens)[i])/M_P;
                    nu_c=calcCyclotronFreq(calcB(el_dens,(hydro_data->temp)[i]));
                    dimlesstheta=calcDimlessTheta( (hydro_data->temp)[i]);
                    //fprintf(fPtr, "B field is: %e at r=%e\n", calcB(el_dens,*(temp+i)), *(r+i));
                    //fflush(fPtr);

                    //printf("Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e\n",*(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta);

                    params[0] = (hydro_data->temp)[i]; //nu_c;
                    params[1]=dimlesstheta;
                    params[2]= el_dens;
                    F.params = &params;
                    
                    //printf("Integrating\n"); //instead integrating from 0 to nu_c
                    status=gsl_integration_qags(&F, 10, nu_c, 0, 1e-2, 10000, w, &ph_dens_calc, &error); //find number of low energy seed photons in the tail of the BB distribution
                    //printf ("error: %s\n", gsl_strerror (status));
                    
                    ph_dens_calc*=hydroElementVolume(hydro_data, i)/(ph_weight_adjusted);
                    
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
            
            //fprintf(fPtr, "dens: %d, photons: %d, adjusted weight: %e\n", *(ph_dens+(j-1)), ph_tot, ph_weight_adjusted);
            //fflush(fPtr);
        }
        
        
        if (block_cnt!=0)
        {
            fprintf(fPtr, "Emitting %d cyclosynchrotron photon(s) with weight %e\n", ph_tot,ph_weight_adjusted );
            fflush(fPtr);
        }
        else
        {
            fprintf(fPtr, "Emitting 0 cyclosynchrotron photons\n" );
            fflush(fPtr);

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
        if (((*ph_orig)[i].weight == 0)) //if photons are null COMPTONIZED_PHOTON photons and not absorbed UNABSORBED_CS_PHOTON photons
        {
            null_ph_count+=1;
        }
    }
    
    
    //allocate memory for that many photons and also allocate memory to hold comoving 4 momentum of each photon and the velocity of the fluid
    //ph_emit=malloc (ph_tot * sizeof (struct photon ));
    p_comv=malloc(4*sizeof(double));
    boost=malloc(4*sizeof(double));
    l_boost=malloc(4*sizeof(double));
    
    if (null_ph_count < ph_tot)
    {
        //if the totoal number of photons to be emitted is larger than the number of null phtons curently in the array, then have to grow the array
        //need to realloc memory to hold the old photon info and the new emitted photon's info
        //fprintf(fPtr, "Emit: Allocating %d space\n", ((*num_ph)+ph_tot-null_ph_count));
        //fflush(fPtr);

        tmp=realloc(*ph_orig, ((*num_ph)+ph_tot-null_ph_count)* sizeof (struct photon )); //may have to look into directly doubling (or *1.5) number of photons each time we need to allocate more memory, can do after looking at profiling for "just enough" memory method
        if (tmp != NULL)
        {
            /* everything ok */
            *ph_orig = tmp;

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
            if (((*ph_orig)[i].weight == 0)   || (i >= *num_ph))
            {
                //preset values for the the newly created spots to hold the emitted phtoons in
                (*ph_orig)[i].weight=0;
                (*ph_orig)[i].nearest_block_index=-1;
                *(null_ph_indexes+j)=i; //save this information so we can use the same syntax for both cases in saving the emitted photon data
                //fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
                //fflush(fPtr);
                j++;
            }
        }
        count_null_indexes=ph_tot; //use this to count the number fo null photons we have actually created, (this can help if we decide to directly double (or *1.5) number of photons each time we need to allocate more memory, then use factor*((*num_ph)+ph_tot)-(*num_ph)
        
        //loop through the original set of photons to see if
        
        //fprintf(fPtr,"Val %d\n", (*(null_ph_indexes+count_null_indexes-1)));
        *num_ph+=net_ph; //update number of photons
        *num_null_ph=ph_tot-null_ph_count; //((*num_ph)+ph_tot)-(*num_ph)-ph_tot; //reserved space - emitted photons-original photons
        //fprintf(fPtr,"New Num PH %d\nNew null hum_ph %d\n", *num_ph, *num_null_ph);
        //fflush(fPtr);
    }
    else
    {
        //otherwise need to find the indexes of these null photons to save the newly emitted photons in them, start searching from the end of the array to efficiently find them
        //dont need to update the number of photons here
        null_ph_indexes=malloc(null_ph_count*sizeof(int));
        j=0;
        for (i=(*num_ph)-1;i>=0;i--)
        {
            if ((*ph_orig)[i].weight == 0)  //if photons are null COMPTONIZED_PHOTON photons and not absorbed UNABSORBED_CS_PHOTON photons
            {
                // if the weight is 0, this is a photons that has been absorbed and is now null
                *(null_ph_indexes+j)=i;
                j++;
                //fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
                //fflush(fPtr);
                
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
        for (i=0;i< hydro_data->num_elements;i++)
        {
            #if DIMENSIONS == 3
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]-0.5*(hydro_data->r2_size)[i]);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]+0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
            #endif

            if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  < rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner < theta_max))
            {
                
                el_dens= ((hydro_data->dens)[i])/M_P;
                nu_c=calcCyclotronFreq(calcB(el_dens,(hydro_data->temp)[i]));
                dimlesstheta=calcDimlessTheta( (hydro_data->temp)[i]);
                max_jnu=2*jnu(nu_c/10, nu_c, dimlesstheta, el_dens);
                
                for(j=0;j<( *(ph_dens+k) ); j++ )
                {
                    //printf("flash_array_idx: %d Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e\n",i, *(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta);

                    fr_dum=nu_c; //set the frequency directly to the cyclotron frequency
                    //fprintf(fPtr, "%lf\n ",fr_dum);
                    //exit(0);
                    #if DIMENSIONS == 2
                        position_phi=gsl_rng_uniform(rand)*2*M_PI;
                    #else
                        position_phi=0;//dont need this in 3D
                    #endif
                    com_v_phi=gsl_rng_uniform(rand)*2*M_PI;
                    com_v_theta=gsl_rng_uniform(rand)*M_PI; //  acos((gsl_rng_uniform(rand)*2)-1) this was for compton scatt, should be isotropic now?
                    //printf("%lf, %lf, %lf\n", position_phi, com_v_phi, com_v_theta);
                    
                    //populate 4 momentum comoving array
                    *(p_comv+0)=PL_CONST*fr_dum/C_LIGHT;
                    *(p_comv+1)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*cos(com_v_phi);
                    *(p_comv+2)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*sin(com_v_phi);
                    *(p_comv+3)=(PL_CONST*fr_dum/C_LIGHT)*cos(com_v_theta);
                    
                    //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
                    #if DIMENSIONS == 3
                        hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], (hydro_data->v2)[i], (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
                    #else
                        //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
                        hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], 0, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
                    #endif
                    (*(boost+0))*=-1;
                    (*(boost+1))*=-1;
                    (*(boost+2))*=-1;
      
                    //printf("%lf, %lf, %lf\n", *(boost+0), *(boost+1), *(boost+2));
                    
                    //boost to lab frame
                    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
                    //printf("Assigning values to struct\n");
                    
                    idx=(*(null_ph_indexes+count_null_indexes-1));
                    //fprintf(fPtr, "Placing photon in index %d\n", idx);
                    (*ph_orig)[idx].p0=(*(l_boost+0));
                    (*ph_orig)[idx].p1=(*(l_boost+1));
                    (*ph_orig)[idx].p2=(*(l_boost+2));
                    (*ph_orig)[idx].p3=(*(l_boost+3));
                    (*ph_orig)[idx].comv_p0=(*(p_comv+0));
                    (*ph_orig)[idx].comv_p1=(*(p_comv+1));
                    (*ph_orig)[idx].comv_p2=(*(p_comv+2));
                    (*ph_orig)[idx].comv_p3=(*(p_comv+3));
                    
                    #if DIMENSIONS == 3
                        hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
                    #else
                        hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
                    #endif
                    (*ph_orig)[idx].r0= cartesian_position_rand_array[0]; //put photons @center of the box with random phi
                    (*ph_orig)[idx].r1= cartesian_position_rand_array[1] ;
                    (*ph_orig)[idx].r2= cartesian_position_rand_array[2]; //y coordinate in flash becomes z coordinate in MCRaT
                    
                    (*ph_orig)[idx].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                    (*ph_orig)[idx].s1=0;
                    (*ph_orig)[idx].s2=0;
                    (*ph_orig)[idx].s3=0;
                    (*ph_orig)[idx].num_scatt=0;
                    (*ph_orig)[idx].weight=ph_weight_adjusted;
                    (*ph_orig)[idx].nearest_block_index=0; //these photons can be scattered
                    (*ph_orig)[idx].type=CS_POOL_PHOTON;
                    //printf("%d\n",ph_tot);
                    ph_tot++; //count how many photons have been emitted
                    count_null_indexes--; //keep track fo the null photon indexes
                    
                    if ((count_null_indexes == 0) || (ph_tot == null_ph_count))
                    {
                        //if count_null_indexes is 0, then all the null photon spaces are filled with emitted photons
                        //if ph_tot is equal to what it used to be
                        i= hydro_data->num_elements;
                        //printf("MCRaT has completed emitting the cyclosynchrotron photons.\n");
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
        
        el_dens= ((hydro_data->dens)[i])/M_P;
        nu_c=calcCyclotronFreq(calcB(el_dens,(hydro_data->temp)[i]));
        
        fr_dum=nu_c; //_scatt; //set the frequency directly to the cyclotron frequency
        //fprintf(fPtr, "%lf %d\n ",fr_dum, (*ph_orig)[scatt_ph_index].nearest_block_index);
        //exit(0);
        #if DIMENSIONS == 2
            position_phi=gsl_rng_uniform(rand)*2*M_PI;
        #else
            position_phi=0;//dont need this in 3D
        #endif
        com_v_phi=gsl_rng_uniform(rand)*2*M_PI;
        com_v_theta=gsl_rng_uniform(rand)*M_PI; //  acos((gsl_rng_uniform(rand)*2)-1) this was for compton scatt, should be isotropic now?
        
        //populate 4 momentum comoving array
        *(p_comv+0)=PL_CONST*fr_dum/C_LIGHT;
        *(p_comv+1)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*cos(com_v_phi);
        *(p_comv+2)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*sin(com_v_phi);
        *(p_comv+3)=(PL_CONST*fr_dum/C_LIGHT)*cos(com_v_theta);
        
        //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
        #if DIMENSIONS == 3
            hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], (hydro_data->v2)[i], (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
        #else
            //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
            hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], 0, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
        #endif
        (*(boost+0))*=-1;
        (*(boost+1))*=-1;
        (*(boost+2))*=-1;
        //printf("%lf, %lf, %lf\n", *(boost+0), *(boost+1), *(boost+2));
        
        //boost to lab frame
        lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
        
        //fprintf(fPtr, "Placing photon in index %d\n", idx);
        (*ph_orig)[idx].p0=(*(l_boost+0));
        (*ph_orig)[idx].p1=(*(l_boost+1));
        (*ph_orig)[idx].p2=(*(l_boost+2));
        (*ph_orig)[idx].p3=(*(l_boost+3));
        (*ph_orig)[idx].comv_p0=(*(p_comv+0));
        (*ph_orig)[idx].comv_p1=(*(p_comv+1));
        (*ph_orig)[idx].comv_p2=(*(p_comv+2));
        (*ph_orig)[idx].comv_p3=(*(p_comv+3));

        #if DIMENSIONS == 3
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
        #else
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
        #endif
        (*ph_orig)[idx].r0= cartesian_position_rand_array[0]; //put photons @center of the box with random phi
        (*ph_orig)[idx].r1= cartesian_position_rand_array[1] ;
        (*ph_orig)[idx].r2= cartesian_position_rand_array[2]; //y coordinate in flash becomes z coordinate in MCRaT

        (*ph_orig)[idx].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
        (*ph_orig)[idx].s1=0;
        (*ph_orig)[idx].s2=0;
        (*ph_orig)[idx].s3=0;
        (*ph_orig)[idx].num_scatt=0;
        (*ph_orig)[idx].weight=(*ph_orig)[scatt_ph_index].weight;
        (*ph_orig)[idx].nearest_block_index=i; //these photons can be scattered
        (*ph_orig)[idx].type=CS_POOL_PHOTON;
        
        //change position of scattered synchrotron photon to be random in the hydro grid
        position_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r0_size)[i])-((hydro_data->r0_size)[i])/2.0; //choose between -size/2 to size/2
        position2_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r1_size)[i])-((hydro_data->r1_size)[i])/2.0;
        #if DIMENSIONS == 3
            position3_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r2_size)[i])-((hydro_data->r2_size)[i])/2.0;
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, (hydro_data->r2)[i]+position3_rand);
        #else
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, position_phi);
        #endif

        //assign random position
        (*ph_orig)[scatt_ph_index].r0=cartesian_position_rand_array[0];
        (*ph_orig)[scatt_ph_index].r1=cartesian_position_rand_array[1];
        (*ph_orig)[scatt_ph_index].r2=cartesian_position_rand_array[2];
        
    }
    //printf("At End of function\n");
    
    {
        free(null_ph_indexes);
    }
    
    //exit(0);
    free(ph_dens); free(p_comv); free(boost); free(l_boost);
    //free(ph_emit);
    
    gsl_integration_workspace_free (w);
    
    return ph_tot;
}

double phAbsCyclosynch(struct photon **ph_orig, int *num_ph, int *num_abs_ph, int *scatt_cyclosynch_num_ph, struct hydro_dataframe *hydro_data, FILE *fPtr) //double *temp, double *dens, FILE *fPtr)
{
    int i=0, count=0, abs_ph_count=0, synch_ph_count=0, num_thread=1;
    int other_count=0;
    #if defined(_OPENMP)
    num_thread=omp_get_num_threads();
    #endif

    double el_dens=0, nu_c=0, abs_count=0;
    //struct photon tmp_ph;//hold temporay photon to move its data
    
    fprintf(fPtr, "In phAbsCyclosynch func begin: abs_ph_count: %d synch_ph_count: %d scatt_cyclosynch_num_ph: %d num_threads: %d\n", abs_ph_count, synch_ph_count, *scatt_cyclosynch_num_ph, num_thread);
    
    *scatt_cyclosynch_num_ph=0;//set thsi equal to 0, to recount in this function and get prepared for the next frame
    
    #pragma omp parallel for num_threads(num_thread) firstprivate(el_dens, nu_c) reduction(+:abs_ph_count)
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight != 0) && ((*ph_orig)[i].nearest_block_index != -1))
        {
            // if the photon isnt a null photon already, see if it should be absorbed
            
            el_dens= (hydro_data->dens)[(*ph_orig)[i].nearest_block_index]/M_P;//(*(dens+(*ph_orig)[i].nearest_block_index))/M_P;
            nu_c=calcCyclotronFreq(calcB(el_dens, (hydro_data->temp)[(*ph_orig)[i].nearest_block_index])); //*(temp+(*ph_orig)[i].nearest_block_index)));
            //printf("photon %d has lab nu %e comv frequency %e and nu_c %e with FLASH grid number %d\n", i, (*ph_orig)[i].p0*C_LIGHT/PL_CONST, (*ph_orig)[i].comv_p0*C_LIGHT/PL_CONST, nu_c, (*ph_orig)[i].nearest_block_index);
            if (((*ph_orig)[i].comv_p0*C_LIGHT/PL_CONST <= nu_c) || ((*ph_orig)[i].type == CS_POOL_PHOTON))
            {
                //if the photon has a frequency less that nu_c, it should be absorbed and becomes a null photon
                //preset values for the the newly created spots to hold the emitted phtoons in;
                
                //if this is a synchrotron photons or photons that have been scattered that were once synch photons in this frame
                //fprintf(fPtr,"photon %d being absorbed\n", i);
                if (((*ph_orig)[i].type != INJECTED_PHOTON) && ((*ph_orig)[i].type != UNABSORBED_CS_PHOTON) )
                {
                    (*ph_orig)[i].weight=0;
                    (*ph_orig)[i].nearest_block_index=-1;
                    abs_ph_count++;
                    
                    if ((*ph_orig)[i].type == CS_POOL_PHOTON)
                    {
                        synch_ph_count++;
                    }
                }
                else
                {
                    //have an injected photon or UNABSORBED_CS_PHOTON (previous COMPTONIZED_PHOTON photon) that has a nu that can be absorbed
                    abs_count+=(*ph_orig)[i].weight;
                    (*ph_orig)[i].p0=-1; //set its energy negative so we know for later analysis that it can't be used and its been absorbed,
                    (*ph_orig)[i].nearest_block_index=-1;
                    //also set the weight equal to 0 since we no longer care about saving it
                    (*ph_orig)[i].weight=0;
                    abs_ph_count++;

                }
            }
            else
            {
                //if the phootn isnt going to be absorbed, see if its a COMPTONIZED_PHOTON photon thats survived and change it to an injected type
                
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
                
                if (((*ph_orig)[i].type == COMPTONIZED_PHOTON) || ((*ph_orig)[i].type == UNABSORBED_CS_PHOTON) )
                {
                    //if the photon is a COMPTONIZED_PHOTON phton (scattered synch photon from the current frame) or a UNABSORBED_CS_PHOTON photon (scattered synch photon) from an old frame
                    //count how many of these there are
                    *scatt_cyclosynch_num_ph+=1;
                }
                
            }
        }
        else
        {
            //see if the photon was a previous INJECTED_PHOTON photon absorbed that we still have to account for in the array
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
                
        //fprintf(fPtr, "photon %d has energy %e and weight %e with FLASH grid number %d\n", i, (*ph_orig)[i].p0*C_LIGHT/1.6e-9, (*ph_orig)[i].weight, (*ph_orig)[i].nearest_block_index);
    }
    //fprintf(fPtr, "In phAbsCyclosynch func: abs_ph_count: %d synch_ph_count: %d scatt_cyclosynch_num_ph: %d\n", abs_ph_count, synch_ph_count, *scatt_cyclosynch_num_ph);
    *num_abs_ph=abs_ph_count; //+synch_ph_count; dont need this
    
    //fprintf(fPtr, "In phAbsCyclosynch func: count before_loop= %d\n", count);

    while (count<*num_ph)
    {
        //overwrite the last few photons to make sure that they are null photons
        (*ph_orig)[count].weight=0;
        (*ph_orig)[count].nearest_block_index=-1;
        //fprintf(fPtr, "photon %d has frequency %e and weight %e with FLASH grid number %d\n", count, (*ph_orig)[count].comv_p0*C_LIGHT/PL_CONST, (*ph_orig)[count].weight, (*ph_orig)[count].nearest_block_index);
        //fflush(fPtr);
        
        count+=1;
    }
    //fprintf(fPtr, "In phAbsCyclosynch func: count after loop= %d\n", count);

    return abs_count;
}


