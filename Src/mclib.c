#include "mcrat.h"

//define constants
const double A_RAD=7.56e-15, C_LIGHT=2.99792458e10, PL_CONST=6.6260755e-27, FINE_STRUCT=7.29735308e-3, CHARGE_EL= 4.8032068e-10;
const double K_B=1.380658e-16, M_P=1.6726231e-24, THOM_X_SECT=6.65246e-25, M_EL=9.1093879e-28 , R_EL=2.817941499892705e-13;

/* Global thread-local RNG pool (declare in header or at file scope) */
static gsl_rng **global_thread_rng = NULL;
static int global_num_threads = 0;

static int compare_pairs(const void *a, const void *b);
static PhotonTimePair *radixSortPairs(PhotonTimePair *pairs, int n, PhotonTimePair *tmp);

void initGlobalThreadRNG(gsl_rng *master_rng, int num_threads)
{
    if (global_thread_rng)
    {
        /* Already initialized, free old ones */
        for (int i = 1; i < global_num_threads; i++)
        {
            gsl_rng_free(global_thread_rng[i]);
        }
        free(global_thread_rng);
    }
    
    global_num_threads = num_threads;
    global_thread_rng = malloc(num_threads * sizeof(gsl_rng *));
    global_thread_rng[0] = master_rng;
    
    const gsl_rng_type *rng_t = gsl_rng_ranlxs0;
    for (int i = 1; i < num_threads; i++)
    {
        global_thread_rng[i] = gsl_rng_alloc(rng_t);
        gsl_rng_set(global_thread_rng[i], gsl_rng_get(master_rng));
    }
}

void freeGlobalThreadRNG(void)
{
    if (global_thread_rng)
    {
        for (int i = 1; i < global_num_threads; i++)
        {
            gsl_rng_free(global_thread_rng[i]);
        }
        free(global_thread_rng);
        global_thread_rng = NULL;
        global_num_threads = 0;
    }
}




void photonInjection(struct photonList *photon_list, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng * rand, FILE *fPtr)
{
    int i=0, block_cnt=0, *ph_dens=NULL, ph_tot=0, j=0,k=0;
    double ph_dens_calc=0.0, fr_dum=1.0, y_dum=0.0, yfr_dum=0.0, fr_max=0, bb_norm=0, position_phi, ph_weight_adjusted, rmin, rmax;
    double com_v_phi, com_v_theta; //comoving phi, theta,
    //double *p_comv=NULL, *boost=NULL, *l_boost=NULL; // comoving 4 momentum for a photon, and boost for photon(to go to lab frame)and pointer to hold array of lorentz boost, to lab frame values //not needed if we are trying to align the memory for these for address contingency with new amd architecture optimization
    double p_comv[4] SIMD_ALIGN;
    double boost[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment
    double l_boost[4] SIMD_ALIGN;
    
    float num_dens_coeff;
    double r_grid_innercorner=0, r_grid_outercorner=0, theta_grid_innercorner=0, theta_grid_outercorner=0;
    double position_rand=0, position2_rand=0, position3_rand=0, cartesian_position_rand_array[3];
    struct photon *ph=NULL, initialized_photon;
    char spect_non_synch='\0';
    
    //first we define spect1 and spec2 depending on the value read in from mc.par
    //need to determine if we need to inject just wien, just bb, just custom, just synch or combo of
    
    
    //find how many blocks are near the injection radius within the angles defined in mc.par, get temperatures and calculate number of photons to allocate memory for 
    //and then rcord which blocks have to have "x" amount of photons injected there
    
    rmin=r_inj - 0.5*C_LIGHT/hydro_data->fps;
    rmax=r_inj + 0.5*C_LIGHT/hydro_data->fps;
    
    //identify if we need to also inject synchrotron photons
    if ((spect == SYNCHROTRON ) || (spect == WIEN_AND_SYNCH ) || (spect == BLACKBODY_AND_SYNCH ) || (spect == CUSTOM_AND_SYNCH ))
    {
        //injecting the synchrotron photons, so pass in 1 for ph_inj_switch
        photonEmitSynch(photon_list, rmin, ph_weight, min_photons, max_photons, 1, theta_min, theta_max, hydro_data, rand, fPtr);
        
        if (spect == WIEN_AND_SYNCH )
        {
            spect_non_synch=WIEN;
        }
        else if (spect == BLACKBODY_AND_SYNCH )
        {
            spect_non_synch=BLACKBODY;
        }
        else if (spect == CUSTOM_AND_SYNCH )
        {
            spect_non_synch=CUSTOM;
        }
        else
        {
            spect_non_synch='\0';
        }
    }
    else
    {
        spect_non_synch = spect;
    }
    
    if (spect_non_synch != '\0')
    {
        //define the number density coeficient, integrate the number density spectrum from 0 to infinity to get this value
        //used to calculate the number density of photons as num_dens_coeff*T_comv^3
        // how should this be defined for the custom spectrum case? -> made it be an input that the user sets
        if (spect_non_synch == WIEN ) //from MCRAT paper, w for wien spectrum
        {
            num_dens_coeff=8.44;
            //printf("in wien spectrum\n");
        }
        else if (spect_non_synch == BLACKBODY)
        {
            num_dens_coeff=20.29; //this is for black body spectrum
            //printf("in BB spectrum");
        }
        else
        {
            //this is for te custom outflow
            num_dens_coeff=PHOTON_DENSITY_COEFF;
        }
        
        
        
        
        for(i=0; i<hydro_data->num_elements; i++)
        {
            #if DIMENSIONS == THREE
                //want inner corner to be close to origin, therfore ned to have abs for 3D cartesian with negative coordinates, shouldnt affect the other geometry systems since theyre all defined from r=0, theta=0, phi=0
                
                //hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]-0.5*(hydro_data->r2_size)[i]);
                //hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]+0.5*(hydro_data->r2_size)[i]);
                
                //therefore do whats below
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, fabs((hydro_data->r0)[i])-0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])-0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])-0.5*(hydro_data->r2_size)[i]);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, fabs((hydro_data->r0)[i])+0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])+0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])+0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
            #endif
            
            //look at all boxes in width delta r=c/fps and within angles we are interested in
            //if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max) && ((hydro_data->r0_size)[i]<1e11) && ((hydro_data->r1_size)[i]<0.09))
            if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max))
            {
                //&& ((hydro_data->r0_size)[i]<1e11) && ((hydro_data->r1)[i]<3.0*3.14/180) is just for testing sph_3d mcrat sim to see if block_cnt is the issue for the 200x normalization issue -> this fixed norm issue, not N_scatt issue when start at frame 0
                // also try injecting photons in frame 1 without above conditions -> didnt fix normalization issue not N_scatt issue
                // also try inj at frame 1 with scale 1e11 -> didnt fixed normalization issue not N_scatt issue
                // also try inj at frame 0 (orig) to see what gets printed for diagnosing CHOMBO refinement levels being an issue
                // try inj at frame 0 with modified if statement and L scale 1e11
                block_cnt++;
                //#if DIMENSIONS == THREE
                //fprintf(fPtr,"rmin %e rmax %e thetamin %e thetamax %e hydro: r0 %e r1 %e r2 %e r0_size %e r1_size %e r2_size %e r_inner %e theta_inner %e r_outer %e theta_outer %e\n", rmin, rmax, theta_min, theta_max, (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i], (hydro_data->r0_size)[i], (hydro_data->r1_size)[i], (hydro_data->r2_size)[i], r_grid_innercorner, theta_grid_innercorner, r_grid_outercorner, theta_grid_outercorner);
                //#else
                //fprintf(fPtr,"rmin %e rmax %e thetamin %e thetamax %e hydro: r0 %e r1 %e r0_size %e r1_size %e r_inner %e theta_inner %e r_outer %e theta_outer %e dens %e\n", rmin, rmax, theta_min, theta_max, (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r0_size)[i], (hydro_data->r1_size)[i], r_grid_innercorner, theta_grid_innercorner, r_grid_outercorner, theta_grid_outercorner, (hydro_data->dens)[i]);
                //#endif
                //fflush(fPtr);
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
            
            for (i=0;i<hydro_data->num_elements;i++)
            {
                //printf("%d\n",i);
                //printf("%e, %e, %e, %e, %e, %e\n", *(r+i),(r_inj - C_LIGHT/fps), (r_inj + C_LIGHT/fps), *(theta+i) , theta_max, theta_min);
                #if DIMENSIONS == THREE
                    //want inner corner to be close to origin, therfore ned to have abs for 3D cartesian with negative coordinates, shouldnt affect the other geometry systems since theyre all defined from r=0, theta=0, phi=0
                    
                    //hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]-0.5*(hydro_data->r2_size)[i]);
                    //hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]+0.5*(hydro_data->r2_size)[i]);
                    
                    //therefore do whats below
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, fabs((hydro_data->r0)[i])-0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])-0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])-0.5*(hydro_data->r2_size)[i]);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, fabs((hydro_data->r0)[i])+0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])+0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])+0.5*(hydro_data->r2_size)[i]);
                    
                #else
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
                #endif
                
                //if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max) && ((hydro_data->r0_size)[i]<1e11) && ((hydro_data->r1_size)[i]<0.09))
                if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max))
                {
                    ph_dens_calc=(4.0/3.0)*hydroElementVolume(hydro_data, i) *(((hydro_data->gamma)[i]*num_dens_coeff*(hydro_data->temp)[i]*(hydro_data->temp)[i]*(hydro_data->temp)[i])/ph_weight_adjusted/hydro_data->fps); //4 comes from L \propto 4p in the limit radiation pressure is greater than the matter energy density and 3 comes from p=u/3, where u is the energy density
                    
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
        ph=malloc (ph_tot * sizeof (struct photon ));
        /*
        p_comv=malloc(4*sizeof(double));
        boost=malloc(3*sizeof(double));
        l_boost=malloc(4*sizeof(double));
        */ //not needed if we are trying to align the memory for these for address contingency with new amd architecture optimization
        
        //go through blocks and assign random energies/locations to proper number of photons
        ph_tot=0;
        k=0;
        //for blackbody injection sampling using Bjorkman and Wood 2001
        double test=0, test_rand1=gsl_rng_uniform_pos(rand), test_rand2=gsl_rng_uniform_pos(rand), test_rand3=gsl_rng_uniform_pos(rand), test_rand4=gsl_rng_uniform_pos(rand), test_rand5=gsl_rng_uniform_pos(rand);
        double test_cnt=0;
        
        for (i=0;i<hydro_data->num_elements;i++)
        {
            #if DIMENSIONS == THREE
                //want inner corner to be close to origin, therfore ned to have abs for 3D cartesian with negative coordinates, shouldnt affect the other geometry systems since theyre all defined from r=0, theta=0, phi=0
                
                //hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]-0.5*(hydro_data->r2_size)[i]);
                //hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], (hydro_data->r2)[i]+0.5*(hydro_data->r2_size)[i]);
                
                //therefore do whats below
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, fabs((hydro_data->r0)[i])-0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])-0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])-0.5*(hydro_data->r2_size)[i]);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, fabs((hydro_data->r0)[i])+0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])+0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])+0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
            #endif
            
            //if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max) && ((hydro_data->r0_size)[i]<1e11) && ((hydro_data->r1_size)[i]<0.09))
            if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  <= rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner <= theta_max))
            {
                
                for(j=0;j<( *(ph_dens+k) ); j++ )
                {
                    //have to get random frequency for the photon comoving frequency
                    if (spect_non_synch == WIEN)
                    {
                        /* old way which also seemed to  be wrong in many ways
                         y_dum=1; //initalize loop
                         yfr_dum=0;
                         while (y_dum>yfr_dum)
                         {
                         fr_dum=gsl_rng_uniform_pos(rand)*6.3e11*((hydro_data->temp)[i]); //in Hz
                         //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
                         y_dum=gsl_rng_uniform_pos(rand);
                         //printf("%lf ",fr_dum);
                         
                         yfr_dum=(1.0/(1.29e31))*pow((fr_dum/((hydro_data->temp)[i])),3.0)/(exp((PL_CONST*fr_dum)/(K_B*((hydro_data->temp)[i]) ))-1); //curve is normalized to maximum
                         }
                         */
                        //now sample from a gamma distribution with a=2 and b=1 (x^2*e^-x) and then convert x=h \nu / k_B*T to frequency \nu
                        //this is due to the fact that we are sampling from the wien photon density spectrum which is the wien spectrum divided by h \nu.
                        //verified with simulations in python to verify this sampled function and transform does actually return a wien function for sufficiently large sample size (of 100000)
                        fr_dum = gsl_ran_gamma(rand, 3.0, 1.0);
                        fr_dum*=K_B*((hydro_data->temp)[i])/PL_CONST;
                        
                    }
                    else if (spect_non_synch == BLACKBODY)
                    {
                        fr_max=(3.31e10)*((hydro_data->temp)[i]);//max frequency of bb photon density spectrum
                        bb_norm=( pow((fr_max),2.0))/gsl_expm1(PL_CONST*fr_max/(K_B*((hydro_data->temp)[i]))); //(exp(PL_CONST*fr_max/(K_B*bb_temp))-1); //find value of bb at fr_max
                        y_dum=1; //initalize loop
                        yfr_dum=0;
                        while (y_dum>yfr_dum)
                        {
                            fr_dum=gsl_rng_uniform_pos(rand)*6.3e11*((hydro_data->temp)[i]); //in Hz
                            //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
                            y_dum=gsl_rng_uniform_pos(rand);
                            
                            yfr_dum=((1.0/bb_norm)* pow((fr_dum),2.0))/gsl_expm1(PL_CONST*fr_dum/(K_B*((hydro_data->temp)[i]))); //(exp(PL_CONST*fr_dum/(K_B*bb_temp))-1); //curve is normalized to vaue of bb @ max frequency
                        }
                        
                    }
                    else
                    {
                        //this is for custom spectrum sampling
                        initialized_photon = custom_photon_sampler(hydro_data, i, rand, fPtr);
                        
                    }
                    //printf("%lf, %lf,%lf,%e \n",(hydro_data->temp)[i],fr_dum, y_dum, yfr_dum);
                    
                    
                    //printf("i: %d freq:%lf\n ",ph_tot, fr_dum);
                    #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
                        position_phi=gsl_rng_uniform(rand)*2*M_PI;
                    #else
                        position_phi=0;//dont need this in 3D
                    #endif
                    com_v_phi=samplePhotonPhi(rand, fPtr); //gsl_rng_uniform(rand)*2*M_PI;
                    //this seemed to produce lab frame spectra with significantly differnet temperatures/shapes than what was expected for wien/blackbody spectra. this is only valid when beta=0, which is limiting case of our anisotropic sampling below
                    //com_v_theta=acos((gsl_rng_uniform(rand)*2)-1);
                    
                    //TODO:what is boost at the start of this loop? it seems undefined
                    //trying to overwrite com_v_theta based on sampling of lab anisotropic angle distribution of photons
                    //see eg Section 3.2.1 @ doi.org/10.1088/0004-637X/807/1/31 & Section 3.5 @ doi.org/10.3847/1538-4357/ac75cb
                    // and section 6.2 here: Nordin Nobuoka, J. 2025, SPIRO: a code that couples Monte Carlo photons to relativistic hydrodynamics - Applications to hot astrophysical plasmas, https://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-368279
                    gsl_vector_view b=gsl_vector_view_array(boost, 3);
                    double beta=gsl_blas_dnrm2(&b.vector);
                    y_dum=1; //initalize loop
                    yfr_dum=0;
                    while (y_dum>yfr_dum)
                    {
                        com_v_theta=2*gsl_rng_uniform_pos(rand)-1; //cos(angle) is from -1 to 1
                        //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
                        y_dum=gsl_rng_uniform_pos(rand);
                        
                        yfr_dum=0.5*(1+beta*com_v_theta); //propability density of angle of photon with respect to fluid motion (doppler boosting factor)
                    }
                    com_v_theta=samplePhotonTheta(boost, rand, fPtr); //acos(com_v_theta);
                    //trying to overwrite com_v_theta based on sampling of lab anisotropic angle distribution of photons
                    
                    
                    //printf("%lf, %lf, %lf\n", position_phi, com_v_phi, com_v_theta);
                    
                    //populate 4 momentum comoving array
                    *(p_comv+0)=PL_CONST*fr_dum/C_LIGHT;
                    *(p_comv+1)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*cos(com_v_phi);
                    *(p_comv+2)=(PL_CONST*fr_dum/C_LIGHT)*sin(com_v_theta)*sin(com_v_phi);
                    *(p_comv+3)=(PL_CONST*fr_dum/C_LIGHT)*cos(com_v_theta);
                    
                    
                    //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
                    #if DIMENSIONS == THREE
                        hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], (hydro_data->v2)[i], (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
                    #elif DIMENSIONS == TWO_POINT_FIVE
                        hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], (hydro_data->v2)[i], (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
                    #else
                        //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
                        hydroVectorToCartesian(boost, (hydro_data->v0)[i], (hydro_data->v1)[i], 0, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
                    #endif
                    (*(boost+0))*=-1;
                    (*(boost+1))*=-1;
                    (*(boost+2))*=-1;
                    
                    //boost to lab frame
                    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
                    //printf("Assignemnt: %e, %e, %e, %e\n", *(l_boost+0), *(l_boost+1), *(l_boost+2),*(l_boost+3));
                    
                    ph[ph_tot].p0=(*(l_boost+0));
                    ph[ph_tot].p1=(*(l_boost+1));
                    ph[ph_tot].p2=(*(l_boost+2));
                    ph[ph_tot].p3=(*(l_boost+3));
                    ph[ph_tot].comv_p0=(*(p_comv+0));
                    ph[ph_tot].comv_p1=(*(p_comv+1));
                    ph[ph_tot].comv_p2=(*(p_comv+2));
                    ph[ph_tot].comv_p3=(*(p_comv+3));
                    
                    //place photons in rand positions within fluid element
                    position_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r0_size)[i])-0.5*((hydro_data->r0_size)[i]); //choose between -size/2 to size/2
                    position2_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r1_size)[i])-0.5*((hydro_data->r1_size)[i]);
                    #if DIMENSIONS == THREE
                        position3_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r2_size)[i])-0.5*((hydro_data->r2_size)[i]);
                        hydroCoordinateToMcratCoordinate(cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, (hydro_data->r2)[i]+position3_rand);
                    #else
                        hydroCoordinateToMcratCoordinate(cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, position_phi);
                    #endif
                    
                    //assign random position
                    ph[ph_tot].r0=cartesian_position_rand_array[0];
                    ph[ph_tot].r1=cartesian_position_rand_array[1];
                    ph[ph_tot].r2=cartesian_position_rand_array[2];
                    
                    //fprintf(fPtr,"%d %e %e %e\n", ph_tot, ph[ph_tot].r0, ph[ph_tot].r1, ph[ph_tot].r2);
                    
                    ph[ph_tot].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                    ph[ph_tot].s1=0;
                    ph[ph_tot].s2=0;
                    ph[ph_tot].s3=0;
                    ph[ph_tot].num_scatt=0;
                    ph[ph_tot].weight=ph_weight_adjusted;
                    ph[ph_tot].nearest_block_index=0;
                    ph[ph_tot].type=INJECTED_PHOTON; //i for injected
                    ph[ph_tot].recalc_properties=1; //set to 1 so we are sure that we calculate tau values later on
                    //printf("%d\n",ph_tot);
                    
                    if ((spect_non_synch != WIEN) && (spect_non_synch != BLACKBODY))
                    {
                        saveUserDefinePhoton((ph+ph_tot), &initialized_photon, hydro_data, i, rand, fPtr);
                    }
                    
                    ph_tot++;
                }
                k++;
            }
        }
        
        //save the whole array to our photon list struct
        //if we have at least 1 photon in the photon_list then we just add the emitted photons, otherwise we have to set the photons in the photon list, this take care of whether synchrotron phtoons have been emitted or not
        if (photon_list->list_capacity >0)
        {
            addToPhotonList(photon_list, ph, ph_tot);
        }
        else
        {
            setPhotonList(photon_list, ph, ph_tot);
        }
    }
    //printf(" %d: %d\n", *(ph_dens+(k-1)), *ph_num);
    free(ph_dens);
    //free(p_comv);free(boost); free(l_boost); //not needed if we are trying to align the memory for these for address contingency with new amd architecture optimization
    free(ph);
    //exit(0);
    
    
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
        gamma=1.0/sqrt(1-beta*beta);
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

int findContainingHydroCell( struct photonList *photon_list, struct hydro_dataframe *hydro_data, int find_nearest_block_switch, gsl_rng * rand, FILE *fPtr)
{
    int i=0, min_index=0, ph_block_index=0, num_thread=1, thread_id=0, num_photons_find_new_element=0;
    bool is_in_block=0; //boolean to determine if the photon is outside of its previously noted block
    double ph_phi=0;
    //double ph_p_comv[4], ph_p[4], fluid_beta[3], photon_hydro_coord[3];
    double ph_p_comv[4] SIMD_ALIGN;
    double ph_p[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment
    double fluid_beta[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment
    double photon_hydro_coord[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment

    struct photon *ph=NULL;

    #if defined(_OPENMP)
        num_thread = omp_get_max_threads();
    #endif
    
    //initialize gsl random number generator fo each thread
    /*
    const gsl_rng_type *rng_t;
    gsl_rng **rng;
    rng_t = gsl_rng_ranlxs0;

    rng = (gsl_rng **) malloc((num_thread ) * sizeof(gsl_rng *));
    rng[0]=rand;

    //#pragma omp parallel for num_threads(nt)
    for(i=1;i<num_thread;i++)
    {
        rng[i] = gsl_rng_alloc (rng_t);
        gsl_rng_set(rng[i],gsl_rng_get(rand));
    }
     */

    //go through each photon and find the blocks around it and then get the distances to all of those blocks and choose the one thats the shortest distance away
    //can optimize here, exchange the for loops and change condition to compare to each of the photons is the radius of the block is .95 (or 1.05) times the min (max) photon radius
    //or just parallelize this part here
    
    #pragma omp parallel for num_threads(num_thread) firstprivate( is_in_block, ph_block_index,  ph_phi, min_index, ph_p_comv, ph_p, fluid_beta, photon_hydro_coord, ph) private(i, thread_id) reduction(+:num_photons_find_new_element)
    for (i=0;i<photon_list->list_capacity; i++)
    {
        ph=getPhoton(photon_list, i);

        //fprintf(fPtr, "%d, %d,%e\n", i, (ph->nearest_block_index), (ph->weight));
        //fflush(fPtr);
        
        if (find_nearest_block_switch==0)
        {
            ph_block_index=ph->nearest_block_index; //if starting a new frame the number of indexes can change and cause a seg fault here
        }
        else
        {
            ph_block_index=0; // therefore if starting a new frame set index=0 to avoid this issue
        }
        
        mcratCoordinateToHydroCoordinate(photon_hydro_coord, ph->r0, ph->r1, ph->r2);//convert the photons coordinate to the hydro sim coordinate system
        
        //printf("ph_x:%e, ph_y:%e\n", ph_x, ph_y);
        
        //if the location of the photon is inside the domain of the hydro simulation then do all of this, otherwise assign huge mfp value so no scattering occurs and the next frame is loaded
        // absorbed photons have ph_block_index=-1, therefore if this value is not less than 0, calulate the mfp properly but doesnt work when go to new frame and find new indexes (will change b/c will get rid of these photons when printing)
        //alternatively make decision based on 0 weight
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
        if (((photon_hydro_coord[1]<(hydro_data->r1_domain)[1]) &&
             (photon_hydro_coord[1]>(hydro_data->r1_domain)[0]) &&
             (photon_hydro_coord[0]<(hydro_data->r0_domain)[1]) &&
             (photon_hydro_coord[0]>(hydro_data->r0_domain)[0])) && (ph->nearest_block_index != -1) ) //can use sorted index to see which photons have been absorbed efficiently before printing and get the indexes
        #else
        if (((photon_hydro_coord[2]<(hydro_data->r2_domain)[1]) &&
             (photon_hydro_coord[2]>(hydro_data->r2_domain)[0]) &&
             (photon_hydro_coord[1]<(hydro_data->r1_domain)[1]) &&
             (photon_hydro_coord[1]>(hydro_data->r1_domain)[0]) &&
             (photon_hydro_coord[0]<(hydro_data->r0_domain)[1]) &&
             (photon_hydro_coord[0]>(hydro_data->r0_domain)[0])) && (ph->nearest_block_index != -1) )
        #endif
        {

            is_in_block=checkInBlock(photon_hydro_coord[0], photon_hydro_coord[1], photon_hydro_coord[2], hydro_data, ph_block_index);
            
            //when rebinning photons can have comoving 4 momenta=0 and nearest_block_index=0 (and block 0 be the actual block the photon is in making it not refind the proper index and reclaulate the comoving 4 momenta) which can make counting synch scattered photons be thrown off, thus take care of this case by forcing the function to recalc things
            #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((ph_block_index==0) && ( (ph->comv_p0)+(ph->comv_p1)+(ph->comv_p2)+(ph->comv_p3) == 0 ) )
                {
                    is_in_block=0; //say that photon is not in the block, force it to recompute things
                }
            #endif
            
            /*
            if (find_nearest_block_switch==0 && is_in_block)
            {
                //keep the saved grid index
                min_index=ph_block_index;
            }
            else
             */
            //if we have purposefully set the find_nearest_block_switch==1 to get the function to recompute the nearest block of the photons
            // or if we have already identified that the photon is not in its currently identified block, we need to find the new hydro cell grid and
            //recompute various quantities
            if (find_nearest_block_switch == 1 || !is_in_block)
            {
                //find the new index of the block closest to the photon
                //min_index=findNearestBlock(array_num,  ph_x,  ph_y,  ph_z,  x,   y,  z); //stop doing this one b/c nearest grid could be one that the photon isnt actually in due to adaptive mesh
            
                //find the new index of the block that the photon is actually in
                min_index=findContainingBlock_grid(photon_hydro_coord[0], photon_hydro_coord[1], photon_hydro_coord[2], hydro_data, fPtr); //(array_num,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy, ph_block_index, find_nearest_block_switch, fPtr);

                ph->nearest_block_index=min_index; //save the index if min_index != -1

                if (min_index != -1)
                {
                    //also recalculate the photons' comoving frequency in this new fluid element
                    ph_p[0]=(ph->p0);
                    ph_p[1]=(ph->p1);
                    ph_p[2]=(ph->p2);
                    ph_p[3]=(ph->p3);
                    
                    #if DIMENSIONS == THREE
                        hydroVectorToCartesian(fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], (hydro_data->v2)[min_index], (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], (hydro_data->r2)[min_index]);
                    #elif DIMENSIONS == TWO_POINT_FIVE
                        ph_phi=atan2((ph->r1), (ph->r0));
                        hydroVectorToCartesian(fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], (hydro_data->v2)[min_index], (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], ph_phi);
                    #else
                        ph_phi=atan2((ph->r1), (ph->r0));
                        //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
                        hydroVectorToCartesian(fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], 0, (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], ph_phi);
                    #endif

                    
                    lorentzBoost(fluid_beta, ph_p, ph_p_comv, 'p', fPtr);
                    
                    (ph->comv_p0)=ph_p_comv[0];
                    (ph->comv_p1)=ph_p_comv[1];
                    (ph->comv_p2)=ph_p_comv[2];
                    (ph->comv_p3)=ph_p_comv[3];

                    #if defined(_OPENMP)
                        thread_id=omp_get_thread_num();
                    #endif

                    //need to also recalculate the optical depth
                    calculateOpticalDepth(ph, hydro_data, global_thread_rng[thread_id], fPtr);
                    if ((ph->recalc_properties)==1)
                    {
                        //if we already needed to recalc the optical depth (due to a scattering or something) else
                        //so just set to 0 since we already got this done
                        (ph->recalc_properties)=0;
                    }


                    num_photons_find_new_element+=1;
                }
                else
                {
                	fprintf(fPtr, "Photon number %d Hydro grid index not found, making sure it doesnt scatter.\n", i);
                }
            
            }
            
        }
        else
        {
            //later on we want to make sure that this photon gets assigned a large mfp so it doesnt erroneously scatter
            //also want to make sure that the absorption optical depth is set to 0
            ph->nearest_block_index=-1;
            //fprintf(fPtr,"Photon %d In ELSE\n", i);
            //exit(0);
            
            
            // Photon is outside the domain or absorbed. Force its time_to_scatter
            // to a large sentinel so it sorts to the back and is never selected
            // by photonEvent within the current frame.
            (ph->time_to_scatter) = 1e12 / C_LIGHT;

        }
        
    }
    
    //free rand number generator
    /*
    for (i=1;i<num_thread;i++)
    {
        gsl_rng_free(rng[i]);
    }
    free(rng);
     */

    //print number of times we had to refind the index of the elemtn photons were located in
    if (find_nearest_block_switch!=0)
    {
        num_photons_find_new_element=0; //force this to be 0 since we forced MCRaT to find the indexes for all the photons here
    }
    
    return num_photons_find_new_element;
    
}

void calcMeanFreePath(struct photonList *photon_list, struct hydro_dataframe *hydro_data, gsl_rng * rand, FILE *fPtr)
{
    int i=0, ph_block_index=0, num_thread=1, thread_id=0;
    double mfp=0, default_mfp=FLT_MAX;
    double rnd_tracker=0;
    struct photon *ph=NULL;
    PhotonTimePair *pairs = photon_list->sort_pairs; //malloc(photon_list->list_capacity * sizeof(PhotonTimePair));
    PhotonTimePair *tmp   = photon_list->sort_tmp; //malloc(photon_list->list_capacity * sizeof(PhotonTimePair));


    #if defined(_OPENMP)
        num_thread = omp_get_max_threads();
    #endif

    //initialize gsl random number generator fo each thread
    /*
    const gsl_rng_type *rng_t;
    gsl_rng **rng;
    rng_t = gsl_rng_ranlxs0;

    rng = (gsl_rng **) malloc((num_thread ) * sizeof(gsl_rng *));
    rng[0]=rand;

    //#pragma omp parallel for num_threads(nt)
    for(i=1;i<num_thread;i++)
    {
        rng[i] = gsl_rng_alloc (rng_t);
        gsl_rng_set(rng[i],gsl_rng_get(rand));
    }
     */

    #pragma omp parallel for num_threads(num_thread) firstprivate(ph_block_index, mfp, rnd_tracker, ph) private(i, thread_id) shared(default_mfp)
    for (i=0;i<photon_list->list_capacity; i++)
    {
        photon_list->sorted_indexes[i]=i; //save  indexes to array to use in qsort

        ph=getPhoton(photon_list, i);
        
        ph_block_index=ph->nearest_block_index;

        //if the location of the photon is inside the domain of the hydro simulation then do all of this, otherwise assign huge mfp value so no scattering occurs and the next frame is loaded
        // absorbed photons have ph_block_index=-1, therefore if this value is not less than 0, calulate the mfp properly but doesnt work when go to new frame and find new indexes (will change b/c will get rid of these photons when printing)
        //alternatively make decision based on 0 weight

        //if ph_block_index!= -1 (know which fluid element photon is in) do all this stuff, otherwise make sure photon doesnt scatter
        if (ph_block_index != -1)
        {
            //fprintf(fPtr,"ph_block Index: %d\n", ph_block_index);


            //put this in to double check that random number is between 0 and 1 (exclusive) because there was a problem with this for parallel case
            rnd_tracker=0;
            #if defined(_OPENMP)
                thread_id=omp_get_thread_num();
            #endif

            if ((ph->recalc_properties)==1)
            {
                //if we need to recalc the optical depth (due to a scattering or something) else then do so
                calculateOpticalDepth(ph, hydro_data, global_thread_rng[thread_id], fPtr);
                (ph->recalc_properties)=0;
            }

            rnd_tracker=gsl_rng_uniform_pos(global_thread_rng[thread_id]);
            //printf("Rnd_tracker: %e Thread number %d \n",rnd_tracker, omp_get_thread_num() );

            //mfp=(-1)*log(rnd_tracker)*(M_P/((n_dens_tmp))/(THOM_X_SECT)); ///(1.0-beta*((n_cosangle)))) ; // the mfp and then multiply it by the ln of a random number to simulate distribution of mean free paths IN COMOV FRAME for reference

            mfp = (-1.0/ph->total_scattering_opacity)*log(rnd_tracker);
        }
        else
        {
            mfp=default_mfp;
            
            #if SYNCHROTRON_SWITCH == ON
                //set this to 0 so the continuous absorption doesnt occur as time steps/path lengths are traversed
                ph->absorption_opacity = 0.0;
            #endif
            
            // Photon is outside the domain or absorbed. Force its time_to_scatter
            // to a large sentinel so it sorts to the back and is never selected
            // by photonEvent within the current frame.
            //this is a bit paranoid, but yolo
            (ph->time_to_scatter) = default_mfp / C_LIGHT;

        }
        
        //save values to use in qsort and to the photon struct itself
        (ph->time_to_scatter)=mfp/C_LIGHT;
        pairs[i].time  = ph->time_to_scatter;
        pairs[i].index = i;
        //fprintf(fPtr,"Photon %d has time %e\n", i, *(all_time_steps+i));
        //fflush(fPtr);

    }
    //exit(0);
    //free rand number generator
    /*
    for (i=1;i<num_thread;i++)
    {
        gsl_rng_free(rng[i]);
    }
    free(rng);
     */

    /* Sort pairs ascending by time — O(N) radix sort */
    PhotonTimePair *sorted_pairs = radixSortPairs(pairs, photon_list->list_capacity, tmp);
    
    /* Extract sorted indices */
    for (i = 0; i < photon_list->list_capacity; i++)
        photon_list->sorted_indexes[i] = sorted_pairs[i].index;

}

static int compare_pairs(const void *a, const void *b)
{
    const PhotonTimePair *pa = (const PhotonTimePair *)a;
    const PhotonTimePair *pb = (const PhotonTimePair *)b;
    return (pa->time > pb->time) - (pa->time < pb->time);
}

/*
 * radixSortPairs
 * --------------
 * Sort an array of PhotonTimePair in ascending order by the time field
 * using an LSD radix sort on the IEEE 754 bit representation of the
 * double time value.
 *
 * Correctness relies on two properties of the time_to_scatter values:
 *   1. Always positive — the sentinel value (FLT_MAX/C_LIGHT) is positive,
 *      and all physical MFP values are positive. For positive doubles,
 *      reinterpreting the bit pattern as uint64_t preserves numerical
 *      order exactly (IEEE 754 guarantee).
 *   2. Always finite — NaN/Inf are excluded by the calcMeanFreePath logic.
 *
 * Algorithm
 * ---------
 * Improvement 1 — Fused histogram:
 *   A single O(N) pass over the input fills all 8 byte-histograms at once,
 *   loading each element's time field from memory exactly once instead of
 *   once per pass (8x reduction in histogram-phase memory traffic).
 *
 * Improvement 2 — Trivial-pass skipping:
 *   After the fused histogram is built, any pass whose histogram has all N
 *   elements in a single bucket is skipped entirely. This is common for the
 *   high bytes of time_to_scatter when all photons share the same dynamical
 *   range (e.g. all times ~ 1e-8 to 1e-5 s). In practice 2-4 of the 8
 *   passes are typically skipped for free.
 *
 * Improvement 3 — Return result pointer:
 *   The function returns a pointer to whichever of the two buffers (pairs
 *   or tmp) holds the sorted result. This is safe regardless of how many
 *   passes were skipped (i.e. regardless of whether the flip count is odd
 *   or even), removing the fragile reliance on exactly 8 passes being even.
 *
 * Parameters
 * ----------
 * pairs : array of PhotonTimePair to sort; may or may not hold the result
 *         on return — use the returned pointer, not this one
 * n     : number of elements
 * tmp   : caller-provided scratch buffer of length >= n
 *
 * Returns
 * -------
 * Pointer to the buffer (either pairs or tmp) that holds the sorted output.
 * The caller must use this pointer, not the original pairs pointer, to read
 * the sorted data.
 */
static PhotonTimePair *radixSortPairs(PhotonTimePair *pairs, int n, PhotonTimePair *tmp)
{
    if (n <= 1) return pairs;

    /* ── Improvement 1: fused histogram ─────────────────────────────────
     * One pass over the input fills all 8 byte-histograms simultaneously.
     * counts[p][b] = number of elements whose byte p equals b.
     * The 8x256 table is 8 KB and fits entirely in L1 cache.            */
    int counts[8][256];
    memset(counts, 0, sizeof(counts));

    for (int i = 0; i < n; i++)
    {
        uint64_t key;
        memcpy(&key, &pairs[i].time, sizeof(uint64_t));

        counts[0][ key        & 0xFF]++;
        counts[1][(key >>  8) & 0xFF]++;
        counts[2][(key >> 16) & 0xFF]++;
        counts[3][(key >> 24) & 0xFF]++;
        counts[4][(key >> 32) & 0xFF]++;
        counts[5][(key >> 40) & 0xFF]++;
        counts[6][(key >> 48) & 0xFF]++;
        counts[7][(key >> 56) & 0xFF]++;
    }

    /* ── Improvements 2 & 3: skip trivial passes, track result buffer ───
     * src always points to the buffer we read from in each scatter pass.
     * dst always points to the buffer we write into.
     * After a real scatter pass the two pointers are swapped.
     * After a trivial (skipped) pass they are NOT swapped, so src stays
     * correct automatically.  The final result is wherever src points.  */
    PhotonTimePair *src = pairs;
    PhotonTimePair *dst = tmp;

    for (int pass = 0; pass < 8; pass++)
    {
        /* ── Improvement 2: check if this pass is trivial ───────────────
         * A pass is trivial when one bucket holds all n elements, meaning
         * every element has the same byte value at this digit position.
         * Sorting them changes nothing, so we skip the scatter entirely.  */
        bool trivial = false;
        for (int b = 0; b < 256; b++)
        {
            if (counts[pass][b] == n) { trivial = true; break; }
        }
        if (trivial) continue;

        /* ── Convert counts to prefix-sum output offsets ─────────────── */
        int offsets[256];
        offsets[0] = 0;
        for (int b = 1; b < 256; b++)
            offsets[b] = offsets[b-1] + counts[pass][b-1];

        /* ── Scatter: read from src, write sorted output to dst ─────── */
        int shift = pass * 8;
        for (int i = 0; i < n; i++)
        {
            uint64_t key;
            memcpy(&key, &src[i].time, sizeof(uint64_t));
            uint8_t byte = (uint8_t)((key >> shift) & 0xFF);
            dst[offsets[byte]++] = src[i];
        }

        /* ── Improvement 3: swap src/dst; result tracks automatically ── */
        PhotonTimePair *swap = src;
        src = dst;
        dst = swap;
    }

    /* src now points to whichever buffer holds the fully sorted result.
     * Return it so the caller reads from the right place regardless of
     * how many passes were skipped.                                      */
    return src;
}


void reverseSortIndexes(void *sorted_indexes, int num_elements, size_t element_size, void *context_array)
{
    /*
    Here, we get the proper call to reverse qsort based on the operating system that we are compiling/running code on
    */
    //printf("before QSORT\n");
    #if (defined _GNU_SOURCE || defined __GNU__ || defined __linux__)
        qsort_r(sorted_indexes, num_elements, element_size,  compare2, context_array);
    #elif (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || defined __FREEBSD__ || defined __BSD__ || defined OpenBSD3_1 || defined OpenBSD3_9)
        qsort_r(sorted_indexes, num_elements, element_size, context_array, compare1);
    #else
        #error Cannot detect operating system
    #endif
}

int compare1(void *ar, const void *a, const void *b)
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

int compare2(const void *a, const void *b, void *ar)
{
    //have 2 compare funcions b/c of changes in qsort_r between BSD and GNU
    //from https://phoxis.org/2012/07/12/get-sorted-index-orderting-of-an-array/
  int aa = *(int *) a;
  int bb = *(int *) b;
  double *arr=NULL;
  arr=ar;
  
return ((arr[aa] > arr[bb]) - (arr[aa] < arr[bb]));
}

int interpolatePropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *n_dens_lab, double *n_vx, double *n_vy, double *n_vz, double *n_temp, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr)
{
    /*
     * THIS FUNCTION IS WRITTEN JUST FOR 2D SIMS AS OF NOW, not used
    */
    /*
    int i=0, j=0, min_index=0, ph_block_index=0, thread_id=0;
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
    double mfp=0,default_mfp=0, beta=0;
        
        
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
    
    default_mfp=1e12;
    #pragma omp parallel for num_threads(num_thread) firstprivate( r, theta,dv, v, all_adjacent_block_indexes, j, left_block_index, right_block_index, top_block_index, bottom_block_index, is_in_block, ph_block_index, ph_x, ph_y, ph_z, ph_phi, min_index, n_dens_lab_tmp,n_vx_tmp, n_vy_tmp, n_vz_tmp, n_temp_tmp, fl_v_x, fl_v_y, fl_v_z, fl_v_norm, ph_v_norm, n_cosangle, mfp, beta, rnd_tracker) private(i) shared(default_mfp )
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
            //min_index=findContainingBlock(array_num,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy, ph_block_index, find_nearest_block_switch, fPtr);
            
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
        #if defined(_OPENMP)
        thread_id=omp_get_thread_num();
        #endif
        
        rnd_tracker=gsl_rng_uniform_pos(rng[thread_id]);
        
        mfp=(-1)*(M_P/((n_dens_lab_tmp))/THOM_X_SECT/(1.0-beta*((n_cosangle))))*log(rnd_tracker) ; //calulate the mfp and then multiply it by the ln of a random number to simulate distribution of mean free paths 
        
        
        #pragma omp critical 
        if ( mfp<default_mfp)
        {
            default_mfp=mfp;
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
            #pragma omp flush(default_mfp)
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
    (*time_step)=default_mfp/C_LIGHT;
    return index;
    */
    return 0;
}


void updatePhotonPosition(struct photonList *photon_list, double t, FILE *fPtr)
{
    //move photons by speed of light
 
    int i=0;
    double old_position=0, new_position=0, divide_p0=0, dl=C_LIGHT*t;
    struct photon *ph=NULL; //pointer to a photon struct
    #if defined(_OPENMP)
        int num_thread=1;
        num_thread = omp_get_max_threads();
    #endif

    
    #pragma omp parallel for num_threads(num_thread) firstprivate(old_position, new_position, divide_p0, ph)
    for (i=0;i<photon_list->list_capacity;i++)
    {
        ph=getPhoton(photon_list, i);
        if ((ph->type != CS_POOL_PHOTON) && (ph->weight != 0))
        {
            old_position= sqrt((ph->r0)*(ph->r0)+(ph->r1)*(ph->r1)+(ph->r2)*(ph->r2)); //uncommented checks since they were not necessary anymore
            
            divide_p0=1.0/(ph->p0);
            
            (ph->r0)+=(ph->p1)*divide_p0*dl; //update x position
            
            (ph->r1)+=(ph->p2)*divide_p0*dl;//update y
            
            (ph->r2)+=(ph->p3)*divide_p0*dl;//update z
            
            new_position= sqrt((ph->r0)*(ph->r0)+(ph->r1)*(ph->r1)+(ph->r2)*(ph->r2));
            /*
            if ((new_position-old_position)/t > C_LIGHT)
            {
                fprintf(fPtr, "PHOTON NUMBER %d IS SUPERLUMINAL. ITS SPEED IS %e c.\n", i, ((new_position-old_position)/t)/C_LIGHT);
            }
            */
            //if ( ph->s0 != 1)
            {
            //	fprintf(fPtr, "PHOTON NUMBER %d DOES NOT HAVE I=1. Instead it is: %e\n", i, ph->s0);
            }
            
            //printf("In update  function: %e, %e, %e, %e, %e, %e, %e\n",(ph->r0), (ph->r1), (ph->r2), t, (ph->p1)/(ph->p0), (ph->p2)/(ph->p0), (ph->p3)/(ph->p0) );
            #if SYNCHROTRON_SWITCH == ON
                applyAbsorption(ph, dl);
            #endif
        }
    }
        
    //printf("In update  function: %e, %e, %e, %e\n",t, ((ph)->p1)/((ph)->p0), ((ph)->p2)/((ph)->p0), ((ph)->p3)/((ph)->p0) );    
    
}






double photonEvent(struct photonList *photon_list, double dt_max, struct hydro_dataframe *hydro_data, int *scattered_ph_index, int *frame_scatt_cnt, int *frame_abs_cnt,  gsl_rng * rand, FILE *fPtr)
{
    //function to perform single photon scattering
    int  i=0, index=0, ph_index=0, event_did_occur=0; //variable event_did_occur is to keep track of wether a scattering or absorption actually occured or not,
    int scattering_subgroup=0; //this is meant for when we have nonthermal electrons to identify which subgroup of electrons we may scatter with
    double scatt_time=0, old_scatt_time=0; //keep track of new time to scatter vs old time to scatter to know how much to incrementally propagate the photons if necessary
    double ph_phi=0, fluid_temp=0;
    /*
    double *ph_p=malloc(4*sizeof(double)); //pointer to hold only photon 4 momentum @ start
    double *el_p_comov=malloc(4*sizeof(double));//pointer to hold the electron 4 momenta in comoving frame
    double *ph_p_comov=malloc(4*sizeof(double));//pointer to hold the comoving photon 4 momenta
    double *fluid_beta=malloc(3*sizeof(double));//pointer to hold fluid velocity vector
    double *negative_fluid_beta=malloc(3*sizeof(double));//pointer to hold negative fluid velocity vector
    double *s=malloc(4*sizeof(double)); //vector to hold the stokes parameters for a given photon
     */ //tried to make this stack allocated
    double ph_p[4] SIMD_ALIGN;
    double el_p_comov[4] SIMD_ALIGN;
    double ph_p_comov[4] SIMD_ALIGN;
    double fluid_beta[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment
    double negative_fluid_beta[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment
    double s[4] SIMD_ALIGN;

    
    struct photon *ph=NULL; //pointer to a photon struct
    bool do_rotation=false; //boolean to help us determine if the stokes parameter needs to be rotated going from lab to fluid frame. We dont need to do this if the fluid is stationary, and if we do then we get a bunch of nans so avoid by setting to false
    
    i=0;
    old_scatt_time=0;
    event_did_occur=0;
    //fprintf(fPtr,"In this function Num_ph %d\n", num_ph);
    //fflush(fPtr);
        
    while (i<photon_list->list_capacity && event_did_occur==0 )
    {
        ph=getPhoton(photon_list, photon_list->sorted_indexes[i]);
        ph_index=photon_list->sorted_indexes[i]; //(*(sorted_indexes+i));
        
        scatt_time= ph->time_to_scatter; //*(all_time_steps+ph_index); //get the time until the photon scatters
        
        //IF THE TIME IS GREATER THAN dt_max dont let the photons positions be updated
        if (scatt_time<dt_max)
        {
            updatePhotonPosition(photon_list, scatt_time-old_scatt_time, fPtr);
        
            //fprintf(fPtr,"i: %d, Photon: %d, Delta t=%e\n", i, ph_index, scatt_time-old_scatt_time);
            //fflush(fPtr);
            
            
            //WHAT IF THE PHOTON MOVES TO A NEW BLOCK BETWEEN WHEN WE CALC MFP AND MOVE IT TO DO THE SCATTERING????
            //it mostly happens at low optical depth, near the photosphere so we would have a large mfp anyways so we probably wouldn't be in this function in that case
            //can also occur for a static fluid with a hard boundary so we account for that here
            // Re-validate position after the move — all photons were advanced
            double photon_hydro_coord[4] SIMD_ALIGN; //this should be 3, but we pad it for memory alignment
            mcratCoordinateToHydroCoordinate(photon_hydro_coord, ph->r0, ph->r1, ph->r2);
            
            #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
                bool ph_still_in_domain =
                    (photon_hydro_coord[0] < (hydro_data->r0_domain)[1]) &&
                    (photon_hydro_coord[0] > (hydro_data->r0_domain)[0]) &&
                    (photon_hydro_coord[1] < (hydro_data->r1_domain)[1]) &&
                    (photon_hydro_coord[1] > (hydro_data->r1_domain)[0]);
            #else
                bool ph_still_in_domain =
                    (photon_hydro_coord[0] < (hydro_data->r0_domain)[1]) &&
                    (photon_hydro_coord[0] > (hydro_data->r0_domain)[0]) &&
                    (photon_hydro_coord[1] < (hydro_data->r1_domain)[1]) &&
                    (photon_hydro_coord[1] > (hydro_data->r1_domain)[0]) &&
                    (photon_hydro_coord[2] < (hydro_data->r2_domain)[1]) &&
                    (photon_hydro_coord[2] > (hydro_data->r2_domain)[0]);
            #endif
            
            //TODO: if we have biasing then we can force scattering in low optical depth regions, so need to make sure we properly capture the properties where the scattering occurs
            //the sorted_indexes gives index of photon with smallest time to potentially scatter then extract the index of the block closest to that photon
            //index can be -1 if the photon has been determined to be absorbed via updatephoton loop above
            index=ph->nearest_block_index;
            
            if (ph_still_in_domain && (index>=0))
            {
                                
                fluid_temp=(hydro_data->temp)[index];
                //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
                
                ph_phi=atan2((ph->r1), ((ph->r0)));
                
                /*
                 if (isnan(ph->r0) || isnan(ph->r1) || isnan(ph->r2))
                 {
                 printf("Not a number\n");
                 }
                 
                 
                 fprintf(fPtr,"ph_phi=%e\n", ph_phi);
                 fflush(fPtr);
                 */
                
                //convert flash coordinated into MCRaT coordinates
                //printf("Getting fluid_beta\n");
                
                #if DIMENSIONS == THREE
                    hydroVectorToCartesian(fluid_beta, (hydro_data->v0)[index], (hydro_data->v1)[index], (hydro_data->v2)[index], (hydro_data->r0)[index], (hydro_data->r1)[index], (hydro_data->r2)[index]);
                #elif DIMENSIONS == TWO_POINT_FIVE
                    hydroVectorToCartesian(fluid_beta, (hydro_data->v0)[index], (hydro_data->v1)[index], (hydro_data->v2)[index], (hydro_data->r0)[index], (hydro_data->r1)[index], ph_phi);
                #else
                    //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
                    hydroVectorToCartesian(fluid_beta, (hydro_data->v0)[index], (hydro_data->v1)[index], 0, (hydro_data->r0)[index], (hydro_data->r1)[index], ph_phi);
                #endif
                
                
                /*
                 fprintf(fPtr,"FLASH v: %e, %e\n", flash_vx,flash_vy);
                 fflush(fPtr);
                 */
                
                //fill in photon 4 momentum
                *(ph_p+0)=(ph->p0);
                *(ph_p+1)=(ph->p1);
                *(ph_p+2)=(ph->p2);
                *(ph_p+3)=(ph->p3);
                
                //first we bring the photon to the fluid's comoving frame
                //already have comoving 4 momentum
                *(ph_p_comov+0)=(ph->comv_p0);
                *(ph_p_comov+1)=(ph->comv_p1);
                *(ph_p_comov+2)=(ph->comv_p2);
                *(ph_p_comov+3)=(ph->comv_p3);
                
                //fill in stokes parameters
                *(s+0)=(ph->s0); //I ==1
                *(s+1)=(ph->s1); //Q/I
                *(s+2)=(ph->s2); //U/I
                *(s+3)=(ph->s3); //V/I
                
                /*
                 if ((ph->type) == COMPTONIZED_PHOTON)
                 {
                 fprintf(fPtr,"Unscattered Photon in Lab frame: %e, %e, %e,%e\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3), (ph->r0), (ph->r1), (ph->r2), *(s+0), *(s+1), *(s+2), *(s+3));
                 fflush(fPtr);
                 fprintf(fPtr,"Fluid Beta: %e, %e, %e\n", *(fluid_beta+0),*(fluid_beta+1), *(fluid_beta+2));
                 fflush(fPtr);
                 }
                 
                 fprintf(fPtr,"Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
                 fflush(fPtr);
                 
                 if ((ph->type) == COMPTONIZED_PHOTON)
                 {
                 fprintf(fPtr, "Before Scattering, In Comov_frame:\n");
                 fflush(fPtr);
                 fprintf(fPtr, "ph_comov: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                 fflush(fPtr);
                 }
                 */
                
                
                //then rotate the stokes plane by some angle such that we are in the stokes coordinat eystsem after the lorentz boost
                #if STOKES_SWITCH == ON
                    //check to see if the fluid is not stationary and we need to do this frame rotation at all, otherwise we get nans
                    do_rotation=(!((*(fluid_beta+0) == 0) && (*(fluid_beta+1) == 0) && (*(fluid_beta+2) == 0)));
                    
                    if (do_rotation)
                    {
                        stokesRotation(fluid_beta, (ph_p+1), (ph_p_comov+1), s, fPtr);
                    }
                #endif
                
                //exit(0);
                //second we generate a thermal/non-thermal electron at the correct temperature
                scattering_subgroup=generateSingleElectron(el_p_comov, fluid_temp, ph_p_comov, ph, rand, fPtr);
                
                /*
                 if ((ph->type) == COMPTONIZED_PHOTON)
                 {
                 fprintf(fPtr,"el_comov: %e, %e, %e,%e\n", *(el_p_comov+0), *(el_p_comov+1), *(el_p_comov+2), *(el_p_comov+3));
                 fflush(fPtr);
                 }
                 */
                
                //third we perform the scattering and save scattered photon 4 monetum in ph_p_comov @ end of function
                event_did_occur=singleScatter(el_p_comov, ph_p_comov, s, rand, fPtr);
                
                /*
                 if ((ph->type) == COMPTONIZED_PHOTON)
                 {
                 fprintf(fPtr,"After Scattering, After Lorentz Boost to Comov frame: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                 fflush(fPtr);
                 }
                 */
                
                if (event_did_occur==1)
                {
                    
                    //we need to make sure that the tau for this photon gets recalculated since we have a new comoving
                    //4 momentum and the photon may no longer be in the same cell (we update the photon position before doing the scattering)
                    //this also means that the photon may be in a completely new cell by the time it scatter though this is unlikely in high density regions
                    //doing this here makes the scattered/unscatterd photons both have to have their properties recalculated
                    (ph->recalc_properties)=1;
                    
                    //fprintf(fPtr,"Within the if!\n");
                    //fflush(fPtr);
                    #if SCATTERING_BIAS_SWITCH == ON
                        // if the scattering bias is 1, we already know that the weight of the nonscattered photon is 0 so can
                        // ignore all of these steps
                        if (ph->scattering_bias[scattering_subgroup] != 1)
                        {
                            double scattered_photon_weight = scatteredPhotonWeight(ph->weight, ph-> scattering_opacity[scattering_subgroup]*(scatt_time-old_scatt_time)*C_LIGHT, ph-> scattering_opacity[scattering_subgroup]*(ph->scattering_bias)[scattering_subgroup] * (scatt_time-old_scatt_time) * C_LIGHT);
                            double unscattered_photon_weight = ph->weight - scattered_photon_weight;
                            
                            if (unscattered_photon_weight<0)
                            {
                                fprintf(fPtr,"The unscattered photon weight is negative!!! \n");
                                fflush(fPtr);
                                exit(1);
                            }

                            
                            //first we set the weight of the scattered photon to be the unscattered weight and then copy it into a new element of the photon_list. This works since none of the fields of the photon struct have been updated based on teh actual scattering yet. That occurs below.
                            ph->weight = unscattered_photon_weight;
                            
                            //add the original to our photon list struct, which does a memcpy into a NULL photon's index
                            //if the photon list has to be expanded, the ph pointer may no longer be valid.
                            //try to get aroudn this by copying the contents of ph pointer to a new photon struct and then pass that in
                            struct photon temp_ph;
                            memcpy(&temp_ph, ph, sizeof(struct photon));
                            addToPhotonList(photon_list, &temp_ph, 1);
                            
                            //now get the address of the scattered photon again incase the photon list was expanded and the original address is no longer valid
                            ph=getPhoton(photon_list, ph_index);
                            
                            //now set the scattered photon weight field  to the correct value
                            ph->weight = scattered_photon_weight;
                            
                            //this will help with SSC not being removed by russian roulette
                            //rebinning should operate on the synch photons and any comptonized photons
                            
                            ph->type = COMPTONIZED_PHOTON;
                            

                        }
                        
                    #endif
                    
                    //if the scattering occured have to uodate the phtoon 4 momentum. if photon didnt scatter nothing changes
                    //fourth we bring the photon back to the lab frame
                    *(negative_fluid_beta+0)=-1*( *(fluid_beta+0));
                    *(negative_fluid_beta+1)=-1*( *(fluid_beta+1));
                    *(negative_fluid_beta+2)=-1*( *(fluid_beta+2));
                    lorentzBoost(negative_fluid_beta, ph_p_comov, ph_p, 'p',  fPtr);
                    
                    /*
                     if ((ph->type) == COMPTONIZED_PHOTON)
                     {
                     fprintf(fPtr,"Scattered Photon in Lab frame: %e, %e, %e,%e\n\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3));
                     fflush(fPtr);
                     }
                     */
                    
                    
                    
                    
                    #if STOKES_SWITCH == ON
                        
                        if (do_rotation)
                        {
                            stokesRotation(negative_fluid_beta, (ph_p_comov+1), (ph_p+1), s, fPtr); //rotate to boost back to lab frame
                        }
                        
                        //save stokes parameters
                        (ph->s0)= *(s+0); //I ==1
                        (ph->s1)= *(s+1);
                        (ph->s2)= *(s+2);
                        (ph->s3)= *(s+3);
                        
                    #endif
                    
                    /*
                    if (((*(ph_p+0))*ENERGY_TO_KEV) > 1e4)
                    {
                        //energy greater than 1e4 keV
                        fprintf(fPtr,"Extremely High Photon Energy!!!!!!!!\n");
                        fflush(fPtr);
                    }
                     */
                    
                    //fprintf(fPtr,"Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
                    //fprintf(fPtr, "Old: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                    
                    
                    //assign the photon its new lab 4 momentum
                    (ph->p0)=(*(ph_p+0));
                    (ph->p1)=(*(ph_p+1));
                    (ph->p2)=(*(ph_p+2));
                    (ph->p3)=(*(ph_p+3));
                    
                    //assign it the comoving frame 4 momentum
                    (ph->comv_p0)=(*(ph_p_comov+0));
                    (ph->comv_p1)=(*(ph_p_comov+1));
                    (ph->comv_p2)=(*(ph_p_comov+2));
                    (ph->comv_p3)=(*(ph_p_comov+3));
                    
                    //printf("Done assigning values to original struct\n");
                    
                    //incremement that photons number of scatterings
                    (ph->num_scatt)+=1;
                    *frame_scatt_cnt+=1; //incrememnt total number of scatterings
                    
                }
            }
            else
            {
                if (index<0)
                {
                    setNullPhoton(photon_list, ph_index);
                }
                else
                {
                    // Photon crossed out of the domain during this position update.
                    // Invalidate its state so it sorts to the back on the next
                    // outer iteration and is not picked for scattering.
                    ph->nearest_block_index = -1;
                    ph->time_to_scatter     = FLT_MAX / C_LIGHT;
                    ph->recalc_properties   = 1;
                }
                
                event_did_occur=1;
           }

        }
        else
        {
            // if the photon scatt_time > dt_max
            //have to adjust the time properly so that the time si now appropriate for the next frame
            scatt_time=dt_max;
            updatePhotonPosition(photon_list, scatt_time-old_scatt_time, fPtr); 
            event_did_occur=1; //set equal to 1 to get out of the loop b/c other subsequent photons will have scatt_time > dt_max
            
        }
    
        old_scatt_time=scatt_time;
        i++;
	}
    //exit(0);
    *scattered_ph_index=ph_index; //save the index of the photon that was scattered
    
    //fprintf(fPtr,"scattered_ph_index: %d %d\n", *scattered_ph_index, (*(sorted_indexes+i-1)));
    //fflush(fPtr);
    
    /*
    free(el_p_comov);
    free(ph_p_comov);
    free(fluid_beta); 
    free(negative_fluid_beta);
    free(ph_p);
    free(s);
    ph_p=NULL;negative_fluid_beta=NULL;ph_p_comov=NULL; el_p_comov=NULL;
    */ //not needed if we create the variables without mallocs
    //retrun total time elapsed to scatter a photon
    return scatt_time;
}

double averagePhotonEnergy(struct photonList *photon_list)
{
    //to calculate weighted photon energy in ergs
    int i=0;
    double e_sum=0, w_sum=0;
    struct photon *ph=NULL;
    
    #pragma omp parallel for reduction(+:e_sum) reduction(+:w_sum)
    for (i=0;i<photon_list->list_capacity;i++)
    {
        ph=getPhoton(photon_list, i);

        #if CYCLOSYNCHROTRON_SWITCH == ON
        if ((ph->weight != 0)) //dont want account for null or absorbed UNABSORBED_CS_PHOTON photons
        #endif
        {
            e_sum+=((ph->p0)*(ph->weight));
            w_sum+=(ph->weight);
        }
    }
    
    return (e_sum*C_LIGHT)/w_sum;
}

void phScattStats(struct photonList *photon_list, int *max, int *min, double *avg, double *r_avg, FILE *fPtr  )
{
    int temp_max=0, temp_min=INT_MAX,  i=0, count=0, count_synch=0, count_comp=0, count_i=0, num_thread=1;
    double sum=0, avg_r_sum=0, avg_r_sum_synch=0, avg_r_sum_comp=0, avg_r_sum_inject=0;
    struct photon *ph=NULL;
    
    #if defined(_OPENMP)
        num_thread = omp_get_max_threads();
    #endif
    
    //printf("Num threads: %d", num_thread);
    #pragma omp parallel for num_threads(num_thread) firstprivate(ph) reduction(min:temp_min) reduction(max:temp_max) reduction(+:sum) reduction(+:avg_r_sum) reduction(+:count)
    for (i=0;i<photon_list->list_capacity;i++)
    {
        ph=getPhoton(photon_list, i);
        
        #if CYCLOSYNCHROTRON_SWITCH == ON
        if ((ph->weight != 0)) //dont want account for null or absorbed UNABSORBED_CS_PHOTON photons
        #endif
        {
            sum+=(ph->num_scatt);
            avg_r_sum+=sqrt((ph->r0)*(ph->r0) + (ph->r1)*(ph->r1) + (ph->r2)*(ph->r2));
            
            //printf("%d %c  %e %e %e %e %e %e\n", i, ph->type, ph->p0, ph->comv_p0, ph->r0, ph->r1, ph->r2, ph->num_scatt);
            
            if ((ph->num_scatt) > temp_max )
            {
                temp_max=(ph->num_scatt);
                //printf("The new max is: %d\n", temp_max);
            }
            
            //if ((i==0) || ((ph->num_scatt)<temp_min))
            if ((ph->num_scatt)<temp_min)
            {
                temp_min=(ph->num_scatt);
                //printf("The new min is: %d\n", temp_min);
            }
            
            if ((ph->type) == INJECTED_PHOTON )
            {
                avg_r_sum_inject+=sqrt((ph->r0)*(ph->r0) + (ph->r1)*(ph->r1) + (ph->r2)*(ph->r2));
                count_i++;
            }
            
            #if CYCLOSYNCHROTRON_SWITCH == ON
            if (((ph->type) == COMPTONIZED_PHOTON) || ((ph->type) == UNABSORBED_CS_PHOTON))
            {
                avg_r_sum_comp+=sqrt((ph->r0)*(ph->r0) + (ph->r1)*(ph->r1) + (ph->r2)*(ph->r2));
                count_comp++;
            }
            #endif
            
            count++;
        }
        
        #if CYCLOSYNCHROTRON_SWITCH == ON
        if ((ph->type) == CS_POOL_PHOTON )
        {
            avg_r_sum_synch+=sqrt((ph->r0)*(ph->r0) + (ph->r1)*(ph->r1) + (ph->r2)*(ph->r2));
            count_synch++;
        }
        #endif
        
    }
    #if CYCLOSYNCHROTRON_SWITCH == ON
        fprintf(fPtr, "In this frame Avg r for i type: %e c and o type: %e and s type: %e\n", avg_r_sum_inject/count_i, avg_r_sum_comp/count_comp, avg_r_sum_synch/count_synch);
    #else
        fprintf(fPtr, "In this frame Avg r for i type: %e \n", avg_r_sum_inject/count_i);
    #endif
    fflush(fPtr);
    //exit(0);
    
    *avg=sum/count;
    *r_avg=avg_r_sum/count;
    *max=temp_max;
    *min=temp_min;
    
}


void phMinMax(struct photonList *photon_list, double *min, double *max, double *min_theta, double *max_theta, FILE *fPtr)
{
    double temp_r_max=0, temp_r_min=DBL_MAX, temp_theta_max=0, temp_theta_min=DBL_MAX;
    int i=0, num_thread=1;
    double ph_r=0, ph_theta=0;
    struct photon *ph=NULL;
    
    #if defined(_OPENMP)
        num_thread = omp_get_max_threads();
    #endif

    
    #pragma omp parallel for num_threads(num_thread) firstprivate(ph_r, ph_theta, ph) reduction(min:temp_r_min) reduction(max:temp_r_max) reduction(min:temp_theta_min) reduction(max:temp_theta_max)
    for (i=0; i<photon_list->list_capacity; i++)
    {
        ph=getPhoton(photon_list, i);
        if (ph->weight != 0)
        {
            ph_r=sqrt((ph->r0)*(ph->r0) + (ph->r1)*(ph->r1) + (ph->r2)*(ph->r2));
            ph_theta=acos((ph->r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            if (ph_r > temp_r_max )
            {
                temp_r_max=ph_r;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, (ph->r0), ph->r1, ph->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_r<temp_r_min)
            {
                temp_r_min=ph_r;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, (ph->r0), ph->r1, ph->r2);
            }
            
            if (ph_theta > temp_theta_max )
            {
                temp_theta_max=ph_theta;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, (ph->r0), ph->r1, ph->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_theta<temp_theta_min)
            {
                temp_theta_min=ph_theta;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, (ph->r0), ph->r1, ph->r2);
            }
        }
    }
    
    *max=temp_r_max;
    *min=temp_r_min;
    *max_theta=temp_theta_max;
    *min_theta=temp_theta_min;
}


void logspace(double start, double stop, int num, double *array)
{
    /**
     * Generate logarithmically-spaced array
     *
     * @param start Starting value (will be 10^start)
     * @param stop Ending value (will be 10^stop)
     * @param num Number of points
     * @param array Output array (must be pre-allocated)
     */

    double step = (stop - start) / (num - 1);
    for (int i = 0; i < num; i++)
    {
        array[i] = pow(10.0, start + i * step);
    }
}

#if SCATTERING_BIAS_SWITCH == ON
    void calculateAverageDimlessTheta(struct hydro_dataframe *hydro_data, FILE *fPtr)
    {
        int i;
        double result=0, temp=0, volume=0;

        for (i=0; i<hydro_data->num_elements; i++)
        {
            temp+=calcDimlessTheta((hydro_data->temp)[i])*hydroElementVolume(hydro_data, i);
            volume+=hydroElementVolume(hydro_data, i);
        }

        hydro_data->average_dimless_theta=temp/volume;
    }
#endif

/*
 * applyRussianRoulette
 * --------------------
 * Cull SYNCH_PHOTON packets whose weights have been reduced below a
 * dynamically computed threshold by the continuous SSA absorption in
 * applyabsorption.
 *
 * Algorithm
 * ---------
 * 1. [PARALLEL] Collect the weights of all live SYNCH_PHOTON packets into
 *    per-thread local arrays, then merge into a single array serially.
 *
 * 2. [SERIAL] Sort the merged weight array and compute the median weight
 *    w_median.  Derive the roulette threshold:
 *
 *      w_thresh = epsilon_rr * w_median
 *
 * 3. [PARALLEL] For each SYNCH_PHOTON packet with w < w_thresh, draw
 *    u ~ Uniform(0,1) using the thread-local RNG from global_thread_rng:
 *
 *      P_surv = w / w_thresh
 *
 *      if u < P_surv:   packet survives, weight boosted to w_thresh
 *      else:            packet killed via setNullPhoton (atomic update
 *                       of photon list counters — see note below)
 *
 *    This is exactly unbiased: E[w_after] = P_surv * w_thresh = w.
 *
 * Thread safety
 * -------------
 * setNullPhoton calls incrementNullPhotonNum which modifies
 * photon_list->num_photons and photon_list->num_null_photons.  These
 * two counter updates are protected by an OpenMP critical section so
 * that concurrent kills from different threads do not produce a race
 * condition.  The photon struct itself is at a unique index in the
 * array, so writing its fields requires no additional locking.
 *
 * Parameters
 * ----------
 * photon_list : the live photon population (modified in-place)
 * epsilon_rr  : dimensionless threshold parameter (e.g. 1e-2);
 *               w_thresh = epsilon_rr * w_median
 * fPtr        : log file
 *
 * Returns
 * -------
 * Number of packets killed (set to NULL_PHOTON).  Returns 0 if there
 * are no SYNCH_PHOTON packets or if no packet falls below w_thresh.
 */
int applyRussianRouletteByType(struct photonList *photon_list,
                         double             epsilon_rr,
                         char               photon_type,
                         FILE              *fPtr)
{
    #if SYNCHROTRON_SWITCH == ON

        int i;
        int n_synch  = 0;
        int n_killed = 0;
        int num_threads = 1;
        int thread_id   = 0;

        #if defined(_OPENMP)
            num_threads = omp_get_max_threads();
        #endif

        if (photon_list == NULL || photon_list->num_photons == 0)
            return 0;

        /* ── Step 1: collect weights of all live SYNCH_PHOTON packets ───────────
         *
         * Each thread fills its own local sub-array to avoid any shared-write
         * contention.  After the parallel region the sub-arrays are merged
         * serially into a single flat array for sorting.
         *
         * We allocate one sub-array of size list_capacity per thread — this
         * is a safe upper bound on the number of SYNCH_PHOTONs any thread
         * could encounter.  The actual count per thread is recorded in
         * thread_n_synch[t].
         */
        int     *thread_n_synch = (int *)    calloc(num_threads, sizeof(int));
        double **thread_weights = (double **)malloc(num_threads * sizeof(double *));

        if (thread_n_synch == NULL || thread_weights == NULL)
        {
            fprintf(fPtr,
                    ">> [applyRussianRoulette] ERROR: malloc failed for "
                    "per-thread arrays. Skipping roulette.\n");
            fflush(fPtr);
            free(thread_n_synch);
            free(thread_weights);
            return 0;
        }

        for (i = 0; i < num_threads; i++)
        {
            thread_weights[i] = (double *)malloc(photon_list->list_capacity
                                                 * sizeof(double));
            if (thread_weights[i] == NULL)
            {
                fprintf(fPtr,
                        ">> [applyRussianRoulette] ERROR: malloc failed for "
                        "thread_weights[%d]. Skipping roulette.\n", i);
                fflush(fPtr);

                /* free already-allocated sub-arrays before returning */
                int k;
                for (k = 0; k < i; k++)
                    free(thread_weights[k]);
                free(thread_weights);
                free(thread_n_synch);
                return 0;
            }
        }

        #pragma omp parallel for                                        \
            num_threads(num_threads)                                    \
            private(i, thread_id)                                       \
            shared(photon_list, thread_weights, thread_n_synch)
        for (i = 0; i < photon_list->list_capacity; i++)
        {
            #if defined(_OPENMP)
                thread_id = omp_get_thread_num();
            #else
                thread_id = 0;
            #endif

            struct photon *ph = getPhoton(photon_list, i);

            //if (ph->type == SYNCH_PHOTON) commented this out since when we have large absorption optical depths, we set the photon weight to be DBL_MIN in applyAbsorption, we definitely want these photons to be removed.
            //We set ph->weight = DBL_MIN instead of 0 also so these photons can potentially be written out instead of causing garbage data being written out in the printPhotons function
            if (ph->type == photon_type && ph->weight > DBL_MIN)
            {
                int local_idx = thread_n_synch[thread_id];
                thread_weights[thread_id][local_idx] = ph->weight;

                #pragma omp atomic
                thread_n_synch[thread_id]++;
            }
        }

        /* Serial merge of per-thread weight sub-arrays into one flat array */
        for (i = 0; i < num_threads; i++)
            n_synch += thread_n_synch[i];

        if (n_synch == 0)
        {
            fprintf(fPtr,
                    ">> [applyRussianRoulette] No %c packets found. "
                    "Skipping roulette.\n", photon_type);
            fflush(fPtr);

            for (i = 0; i < num_threads; i++)
                free(thread_weights[i]);
            free(thread_weights);
            free(thread_n_synch);
            return 0;
        }

        double *weights = (double *)malloc(n_synch * sizeof(double));
        if (weights == NULL)
        {
            fprintf(fPtr,
                    ">> [applyRussianRoulette] ERROR: malloc failed for merged "
                    "weights array. Skipping roulette.\n");
            fflush(fPtr);

            for (i = 0; i < num_threads; i++)
                free(thread_weights[i]);
            free(thread_weights);
            free(thread_n_synch);
            return 0;
        }

        int offset = 0;
        for (i = 0; i < num_threads; i++)
        {
            memcpy(weights + offset, thread_weights[i],
                   thread_n_synch[i] * sizeof(double));
            offset += thread_n_synch[i];
            free(thread_weights[i]);
        }
        free(thread_weights);
        free(thread_n_synch);

        /* ── Step 2: sort and compute median weight (serial) ─────────────────────
         *
         * gsl_sort sorts in-place in ascending order.
         * For odd  n_synch: median = weights[n_synch/2]
         * For even n_synch: median = mean of two central values
         */
        gsl_sort(weights, 1, (size_t)n_synch);

        double w_median;
        if (n_synch % 2 == 1)
            w_median = weights[n_synch / 2];
        else
            w_median = 0.5 * (weights[n_synch / 2 - 1] + weights[n_synch / 2]);

        free(weights);

        if (w_median <= 0.0)
        {
            fprintf(fPtr,
                    ">> [applyRussianRoulette] WARNING: median %c "
                    "weight is <= 0 (w_median = %.3e). Skipping roulette.\n",
                    photon_type, w_median);
            fflush(fPtr);
            return 0;
        }

        double w_thresh = epsilon_rr * w_median;

        fprintf(fPtr,
                ">> [applyRussianRoulette] n_synch=%d  w_median=%.3e  "
                "epsilon_rr=%.3e  w_thresh=%.3e\n",
                n_synch, w_median, epsilon_rr, w_thresh);
        fflush(fPtr);

        /* ── Step 3: roulette pass (parallel) ────────────────────────────────────
         *
         * Each photon slot is owned by exactly one thread, so reading and
         * writing ph->weight and ph->type is race-free without locking.
         *
         * The only shared-state modification is the photon list counters
         * (num_photons, num_null_photons) inside setNullPhoton via
         * incrementNullPhotonNum.  These are protected by a critical section.
         *
         * n_killed is accumulated with an atomic update.
         */
        #pragma omp parallel for                                        \
            num_threads(num_threads)                                    \
            private(i, thread_id)                                       \
            shared(photon_list, w_thresh, n_killed)
        for (i = 0; i < photon_list->list_capacity; i++)
        {
            #if defined(_OPENMP)
                thread_id = omp_get_thread_num();
            #else
                thread_id = 0;
            #endif

            struct photon *ph = getPhoton(photon_list, i);

            if (ph->type == photon_type && ph->weight < w_thresh)
            {
                double p_surv = ph->weight / w_thresh;
                double u      = gsl_rng_uniform_pos(global_thread_rng[thread_id]);

                if (u < p_surv)
                {
                    /* Packet survives: boost weight to threshold.
                     * Only this thread touches this slot — no lock needed. */
                    ph->weight = w_thresh;
                }
                else
                {
                    /* Packet killed.
                     * setNullPhoton writes to this slot (race-free) but also
                     * calls incrementNullPhotonNum which modifies the shared
                     * list counters — protect with a critical section.       */
                    #pragma omp critical (roulette_kill)
                    {
                        setNullPhoton(photon_list, i);
                    }

                    #pragma omp atomic
                    n_killed++;
                }
            }
        }

        fprintf(fPtr,
                ">> [applyRussianRoulette] Complete: %d of %d %c "
                "packets killed (%.1f%%).\n\n",
                n_killed, n_synch, photon_type,
                (n_synch > 0) ? 100.0 * n_killed / n_synch : 0.0);
        fflush(fPtr);

        return n_killed;

    #else
        (void)photon_list;
        (void)epsilon_rr;
        (void)fPtr;
        return 0;
    #endif
}

int applyRussianRoulette(struct photonList *photon_list,
                         double             epsilon_rr,
                         FILE              *fPtr)
{
    #if SYNCHROTRON_SWITCH == ON
        int n_killed=0;

        n_killed = applyRussianRouletteByType(photon_list, epsilon_rr, SYNCH_PHOTON, fPtr);
    
        //injected photons can scatter with power law electrons and get very small weights
        n_killed = applyRussianRouletteByType(photon_list, epsilon_rr, INJECTED_PHOTON, fPtr);

        //want to get rid of the SSC photons that have extremely low weights b/c the synch photons with low weights got scattered before they were removed. If we use the same epsilon_rr then it will be hard to capture the high energy portion of the SSC spectrum therefore reduce the epsilon_rr by some other factor of 1e8. This number was tested based on the broken powerlaw SSC test in the RAIKOU paper (Figure 11 of Kawashima et al 2023 ApJ 949 101)
        //also put this last so it actually returns the number of comptonized photons that have been absorbed so we can know when to rebin
        n_killed = applyRussianRouletteByType(photon_list, 1e-8*epsilon_rr, COMPTONIZED_PHOTON, fPtr);

        return n_killed;

    #else
        (void)photon_list;
        (void)epsilon_rr;
        (void)fPtr;
        return 0;
    #endif
}
