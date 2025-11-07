#include "mcrat.h"


//define constants
const double A_RAD=7.56e-15, C_LIGHT=2.99792458e10, PL_CONST=6.6260755e-27, FINE_STRUCT=7.29735308e-3, CHARGE_EL= 4.8032068e-10;
const double K_B=1.380658e-16, M_P=1.6726231e-24, THOM_X_SECT=6.65246e-25, M_EL=9.1093879e-28 , R_EL=2.817941499892705e-13;


void photonInjection(struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng * rand, FILE *fPtr)
{
    int i=0, block_cnt=0, *ph_dens=NULL, ph_tot=0, j=0,k=0;
    double ph_dens_calc=0.0, fr_dum=0.0, y_dum=0.0, yfr_dum=0.0, fr_max=0, bb_norm=0, position_phi, ph_weight_adjusted, rmin, rmax;
    double com_v_phi, com_v_theta, *p_comv=NULL, *boost=NULL; //comoving phi, theta, comoving 4 momentum for a photon, and boost for photon(to go to lab frame)
    double *l_boost=NULL; //pointer to hold array of lorentz boost, to lab frame, values
    float num_dens_coeff;
    double r_grid_innercorner=0, r_grid_outercorner=0, theta_grid_innercorner=0, theta_grid_outercorner=0;
    double position_rand=0, position2_rand=0, position3_rand=0, cartesian_position_rand_array[3];
    
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
    
    rmin=r_inj - 0.5*C_LIGHT/hydro_data->fps;
    rmax=r_inj + 0.5*C_LIGHT/hydro_data->fps;
    
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
                ph_dens_calc=(4.0/3.0)*hydroElementVolume(hydro_data, i) *(((hydro_data->gamma)[i]*num_dens_coeff*(hydro_data->temp)[i]*(hydro_data->temp)[i]*(hydro_data->temp)[i])/ph_weight_adjusted); //4 comes from L \propto 4p in the limit radiation pressure is greater than the matter energy density and 3 comes from p=u/3, where u is the energy density
                
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
                if (spect=='w')
                {
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
                }
                else
                {
                        /* old way
                        fr_max=(5.88e10)*((hydro_data->temp)[i]);//(C_LIGHT*(*(temps+i)))/(0.29); //max frequency of bb
                        bb_norm=(PL_CONST*fr_max * pow((fr_max/C_LIGHT),2.0))/(exp(PL_CONST*fr_max/(K_B*((hydro_data->temp)[i])))-1); //find value of bb at fr_max
                        yfr_dum=((1.0/bb_norm)*PL_CONST*fr_dum * pow((fr_dum/C_LIGHT),2.0))/(exp(PL_CONST*fr_dum/(K_B*((hydro_data->temp)[i])))-1); //curve is normalized to vaue of bb @ max frequency
                        */
                        
                        test=0;
                        test_rand1=gsl_rng_uniform_pos(rand);
                        test_rand2=gsl_rng_uniform_pos(rand);
                        test_rand3=gsl_rng_uniform_pos(rand);
                        test_rand4=gsl_rng_uniform_pos(rand);
                        test_rand5=gsl_rng_uniform_pos(rand);
                        test_cnt=0;
                        while (test<M_PI*M_PI*M_PI*M_PI*test_rand1/90.0)
                        {
                            test_cnt+=1;
                            test+=1/(test_cnt*test_cnt*test_cnt*test_cnt);
                        }
                        fr_dum=-log(test_rand2*test_rand3*test_rand4*test_rand5)/test_cnt;
                        fr_dum*=K_B*((hydro_data->temp)[i])/PL_CONST;
                        y_dum=0; yfr_dum=1;
                        
                }
                    //printf("%lf, %lf,%lf,%e \n",(*(temps+i)),fr_dum, y_dum, yfr_dum);
                    
                
                //printf("i: %d freq:%lf\n ",ph_tot, fr_dum);
                #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
                    position_phi=gsl_rng_uniform(rand)*2*M_PI;
                #else
                    position_phi=0;//dont need this in 3D
                #endif
               com_v_phi=gsl_rng_uniform(rand)*2*M_PI;
               com_v_theta=acos((gsl_rng_uniform(rand)*2)-1);
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
               
                (*ph)[ph_tot].p0=(*(l_boost+0));
                (*ph)[ph_tot].p1=(*(l_boost+1));
                (*ph)[ph_tot].p2=(*(l_boost+2));
                (*ph)[ph_tot].p3=(*(l_boost+3));
                (*ph)[ph_tot].comv_p0=(*(p_comv+0));
                (*ph)[ph_tot].comv_p1=(*(p_comv+1));
                (*ph)[ph_tot].comv_p2=(*(p_comv+2));
                (*ph)[ph_tot].comv_p3=(*(p_comv+3));
                
                //place photons in rand positions within fluid element
                position_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r0_size)[i])-0.5*((hydro_data->r0_size)[i]); //choose between -size/2 to size/2
                position2_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r1_size)[i])-0.5*((hydro_data->r1_size)[i]);
                #if DIMENSIONS == THREE
                    position3_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r2_size)[i])-0.5*((hydro_data->r2_size)[i]);
                    hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, (hydro_data->r2)[i]+position3_rand);
                #else
                    hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, position_phi);
                #endif
                
                //assign random position
                (*ph)[ph_tot].r0=cartesian_position_rand_array[0];
                (*ph)[ph_tot].r1=cartesian_position_rand_array[1];
                (*ph)[ph_tot].r2=cartesian_position_rand_array[2];
                
                //fprintf(fPtr,"%d %e %e %e\n", ph_tot, (*ph)[ph_tot].r0, (*ph)[ph_tot].r1, (*ph)[ph_tot].r2);
                
                (*ph)[ph_tot].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1 
                (*ph)[ph_tot].s1=0;
                (*ph)[ph_tot].s2=0;
                (*ph)[ph_tot].s3=0;
                (*ph)[ph_tot].num_scatt=0;
                (*ph)[ph_tot].weight=ph_weight_adjusted;
                (*ph)[ph_tot].nearest_block_index=0;
                (*ph)[ph_tot].type=INJECTED_PHOTON; //i for injected
                //printf("%d\n",ph_tot);
                ph_tot++;
            }
            k++;
        }
    }
    
    *ph_num=ph_tot; //save number of photons
    //printf(" %d: %d\n", *(ph_dens+(k-1)), *ph_num);
    free(ph_dens); free(p_comv);free(boost); free(l_boost);
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

int findNearestPropertiesAndMinMFP( struct photon *ph, int num_ph, double *all_time_steps, int *sorted_indexes, struct hydro_dataframe *hydro_data, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr)
{
    int i=0, min_index=0, ph_block_index=0, num_thread=1, thread_id=0;
    double ph_x=0, ph_y=0, ph_phi=0, ph_z=0, ph_r=0, ph_theta=0;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates

    double ph_v_norm=0, fl_v_norm=0, synch_x_sect=0;
    double n_cosangle=0, n_dens_lab_tmp=0,n_vx_tmp=0, n_vy_tmp=0, n_vz_tmp=0, n_temp_tmp=0 ;
    double rnd_tracker=0, n_dens_min=0, n_vx_min=0, n_vy_min=0, n_vz_min=0, n_temp_min=0;
    #if defined(_OPENMP)
    num_thread=omp_get_num_threads(); //default is one above if theres no openmp usage
    #endif
    bool is_in_block=0; //boolean to determine if the photon is outside of its previously noted block
    
    int index=0, num_photons_find_new_element=0;
    double mfp=0, default_mfp=0, beta=0;
    double el_p[4];
    double ph_p_comv[4], ph_p[4], fluid_beta[3], photon_hydro_coord[3];

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
    #pragma omp parallel for num_threads(num_thread) firstprivate( is_in_block, ph_block_index, ph_x, ph_y, ph_z, ph_phi, ph_r, min_index, n_dens_lab_tmp,n_vx_tmp, n_vy_tmp, n_vz_tmp, n_temp_tmp, fl_v_x, fl_v_y, fl_v_z, fl_v_norm, ph_v_norm, n_cosangle, mfp, beta, rnd_tracker, ph_p_comv, el_p, ph_p, fluid_beta) private(i) shared(default_mfp ) reduction(+:num_photons_find_new_element)
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
        
        mcratCoordinateToHydroCoordinate(&photon_hydro_coord, (ph+i)->r0, (ph+i)->r1, (ph+i)->r2);//convert the photons coordinate to the hydro sim coordinate system
        
        //printf("ph_x:%e, ph_y:%e\n", ph_x, ph_y);
        
        //if the location of the photon is inside the domain of the hydro simulation then do all of this, otherwise assign huge mfp value so no scattering occurs and the next frame is loaded
        // absorbed photons have ph_block_index=-1, therefore if this value is not less than 0, calulate the mfp properly but doesnt work when go to new frame and find new indexes (will change b/c will get rid of these photons when printing)
        //alternatively make decision based on 0 weight
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
        if (((photon_hydro_coord[1]<(hydro_data->r1_domain)[1]) &&
             (photon_hydro_coord[1]>(hydro_data->r1_domain)[0]) &&
             (photon_hydro_coord[0]<(hydro_data->r0_domain)[1]) &&
             (photon_hydro_coord[0]>(hydro_data->r0_domain)[0])) && ((ph+i)->nearest_block_index != -1) ) //can use sorted index to see which photons have been absorbed efficiently before printing and get the indexes
        #else
        if (((photon_hydro_coord[2]<(hydro_data->r2_domain)[1]) &&
             (photon_hydro_coord[2]>(hydro_data->r2_domain)[0]) &&
             (photon_hydro_coord[1]<(hydro_data->r1_domain)[1]) &&
             (photon_hydro_coord[1]>(hydro_data->r1_domain)[0]) &&
             (photon_hydro_coord[0]<(hydro_data->r0_domain)[1]) &&
             (photon_hydro_coord[0]>(hydro_data->r0_domain)[0])) && ((ph+i)->nearest_block_index != -1) )
        #endif
        {

            is_in_block=checkInBlock(photon_hydro_coord[0], photon_hydro_coord[1], photon_hydro_coord[2], hydro_data, ph_block_index);
            
            //when rebinning photons can have comoving 4 momenta=0 and nearest_block_index=0 (and block 0 be the actual block the photon is in making it not refind the proper index and reclaulate the comoving 4 momenta) which can make counting synch scattered photons be thrown off, thus take care of this case by forcing the function to recalc things
            #if CYCLOSYNCHROTRON_SWITCH == ON
                if ((ph_block_index==0) && ( ((ph+i)->comv_p0)+((ph+i)->comv_p1)+((ph+i)->comv_p2)+((ph+i)->comv_p3) == 0 ) )
                {
                    is_in_block=0; //say that photon is not in the block, force it to recompute things
                }
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
                min_index=findContainingBlock(photon_hydro_coord[0], photon_hydro_coord[1], photon_hydro_coord[2], hydro_data, fPtr); //(array_num,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy, ph_block_index, find_nearest_block_switch, fPtr);
                
                if (min_index != -1)
                {
                    (ph+i)->nearest_block_index=min_index; //save the index if min_index != -1
                    
                    //also recalculate the photons' comoving frequency in this new fluid element
                    ph_p[0]=((ph+i)->p0);
                    ph_p[1]=((ph+i)->p1);
                    ph_p[2]=((ph+i)->p2);
                    ph_p[3]=((ph+i)->p3);
                    
                    #if DIMENSIONS == THREE
                        hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], (hydro_data->v2)[min_index], (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], (hydro_data->r2)[min_index]);
                    #elif DIMENSIONS == TWO_POINT_FIVE
                        ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
                        hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], (hydro_data->v2)[min_index], (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], ph_phi);
                    #else
                        ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
                        //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
                        hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], 0, (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], ph_phi);
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
                (n_dens_lab_tmp)= (hydro_data->dens_lab)[min_index];//(*(dens_lab+min_index));
                (n_temp_tmp)= (hydro_data->temp)[min_index];//(*(temp+min_index));
                
                #if DIMENSIONS == THREE
                    hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], (hydro_data->v2)[min_index], (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], (hydro_data->r2)[min_index]);
                #elif DIMENSIONS == TWO_POINT_FIVE
                    ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
                    hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], (hydro_data->v2)[min_index], (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], ph_phi);
                #else
                    ph_phi=atan2(((ph+i)->r1), ((ph+i)->r0));
                    //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
                    hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[min_index], (hydro_data->v1)[min_index], 0, (hydro_data->r0)[min_index], (hydro_data->r1)[min_index], ph_phi);
                #endif
                
                fl_v_x=fluid_beta[0];
                fl_v_y=fluid_beta[1];
                fl_v_z=fluid_beta[2];
                
                fl_v_norm=sqrt(fl_v_x*fl_v_x+fl_v_y*fl_v_y+fl_v_z*fl_v_z);
                ph_v_norm=sqrt(((ph+i)->p1)*((ph+i)->p1)+((ph+i)->p2)*((ph+i)->p2)+((ph+i)->p3)*((ph+i)->p3));
        
                //(*(n_cosangle+i))=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //find cosine of the angle between the photon and the fluid velocities via a dot product
                n_cosangle=((fl_v_x* ((ph+i)->p1))+(fl_v_y* ((ph+i)->p2))+(fl_v_z* ((ph+i)->p3)))/(fl_v_norm*ph_v_norm ); //make 1 for cylindrical otherwise its undefined
                
                beta=sqrt(1.0-1.0/((hydro_data->gamma)[min_index]*(hydro_data->gamma)[min_index]));
                
                //put this in to double check that random number is between 0 and 1 (exclusive) because there was a problem with this for parallel case
                rnd_tracker=0;
                #if defined(_OPENMP)
                thread_id=omp_get_thread_num();
                #endif
                
                rnd_tracker=gsl_rng_uniform_pos(rng[thread_id]);
                //printf("Rnd_tracker: %e Thread number %d \n",rnd_tracker, omp_get_thread_num() );
        
                //mfp=(-1)*log(rnd_tracker)*(M_P/((n_dens_tmp))/(THOM_X_SECT)); ///(1.0-beta*((n_cosangle)))) ; // the mfp and then multiply it by the ln of a random number to simulate distribution of mean free paths IN COMOV FRAME for reference
                mfp=(-1)*(M_P/((n_dens_lab_tmp))/THOM_X_SECT/(1.0-beta*n_cosangle))*log(rnd_tracker) ;
                
                
            }
            else
            {
                mfp=default_mfp;
            }
        }
         else
        {
            mfp=default_mfp;
            //fprintf(fPtr,"Photon %d In ELSE\n", i);
            //exit(0);
        }
        
        *(all_time_steps+i)=mfp/C_LIGHT;
        //fprintf(fPtr,"Photon %d has time %e\n", i, *(all_time_steps+i));
        //fflush(fPtr);
        
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
    
    //print number of times we had to refind the index of the elemtn photons were located in
    if (find_nearest_block_switch!=0)
    {
        num_photons_find_new_element=0; //force this to be 0 since we forced MCRaT to find the indexes for all the photons here
    }
    
    return num_photons_find_new_element;
    
}

int compare1 (void *ar, const void *a, const void *b)
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


void updatePhotonPosition(struct photon *ph, int num_ph, double t, FILE *fPtr)
{
    //move photons by speed of light
 
    int i=0;
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    double old_position=0, new_position=0, divide_p0=0;
    
    
    #pragma omp parallel for num_threads(num_thread) firstprivate(old_position, new_position, divide_p0)
    for (i=0;i<num_ph;i++)
    {
        if (((ph+i)->type != CS_POOL_PHOTON) && ((ph+i)->weight != 0))
        {
            old_position= sqrt(((ph+i)->r0)*((ph+i)->r0)+((ph+i)->r1)*((ph+i)->r1)+((ph+i)->r2)*((ph+i)->r2)); //uncommented checks since they were not necessary anymore
            
            divide_p0=1.0/((ph+i)->p0);
            
            ((ph+i)->r0)+=((ph+i)->p1)*divide_p0*C_LIGHT*t; //update x position
            
            ((ph+i)->r1)+=((ph+i)->p2)*divide_p0*C_LIGHT*t;//update y
            
            ((ph+i)->r2)+=((ph+i)->p3)*divide_p0*C_LIGHT*t;//update z
            
            new_position= sqrt(((ph+i)->r0)*((ph+i)->r0)+((ph+i)->r1)*((ph+i)->r1)+((ph+i)->r2)*((ph+i)->r2));
            /*
            if ((new_position-old_position)/t > C_LIGHT)
            {
                fprintf(fPtr, "PHOTON NUMBER %d IS SUPERLUMINAL. ITS SPEED IS %e c.\n", i, ((new_position-old_position)/t)/C_LIGHT);
            }
            */
            //if ( (ph+i)->s0 != 1)
            {
            //	fprintf(fPtr, "PHOTON NUMBER %d DOES NOT HAVE I=1. Instead it is: %e\n", i, (ph+i)->s0);
            }
            
            //printf("In update  function: %e, %e, %e, %e, %e, %e, %e\n",((ph+i)->r0), ((ph+i)->r1), ((ph+i)->r2), t, ((ph+i)->p1)/((ph+i)->p0), ((ph+i)->p2)/((ph+i)->p0), ((ph+i)->p3)/((ph+i)->p0) );
        }
    }
        
    //printf("In update  function: %e, %e, %e, %e\n",t, ((ph)->p1)/((ph)->p0), ((ph)->p2)/((ph)->p0), ((ph)->p3)/((ph)->p0) );    
    
}






double photonEvent(struct photon *ph, int num_ph, double dt_max, double *all_time_steps, int *sorted_indexes, struct hydro_dataframe *hydro_data, int *scattered_ph_index, int *frame_scatt_cnt, int *frame_abs_cnt,  gsl_rng * rand, FILE *fPtr)//(struct photon *ph, int num_ph, double dt_max, double *all_time_steps, int *sorted_indexes, double *all_flash_vx, double *all_flash_vy, double *all_flash_vz, double *all_fluid_temp, int *scattered_ph_index, int *frame_scatt_cnt, int *frame_abs_cnt, gsl_rng * rand, FILE *fPtr)
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
            
            
            //WHAT IF THE PHOTON MOVES TO A NEW BLOCK BETWEEN WHEN WE CALC MFP AND MOVE IT TO DO THE SCATTERING????
            //it mostly happens at low optical depth, near the photosphere so we would have a large mfp anyways so we probably wouldn't be in this function in that case
            index=(ph+ph_index)->nearest_block_index; //the sorted_indexes gives index of photon with smallest time to potentially scatter then extract the index of the block closest to that photon
    
            fluid_temp=(hydro_data->temp)[index];
            //if (strcmp(DIM_SWITCH, dim_3d_str)==0)
    
            ph_phi=atan2(((ph+ph_index)->r1), (((ph+ph_index)->r0)));
            
            /*
            if (isnan((ph+ph_index)->r0) || isnan((ph+ph_index)->r1) || isnan((ph+ph_index)->r2))
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
            *(ph_p+0)=((ph+ph_index)->p0);
            *(ph_p+1)=((ph+ph_index)->p1);
            *(ph_p+2)=((ph+ph_index)->p2);
            *(ph_p+3)=((ph+ph_index)->p3);
            
            //first we bring the photon to the fluid's comoving frame
            //already have comoving 4 momentum
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
            if (((ph+ph_index)->type) == COMPTONIZED_PHOTON)
            {
            fprintf(fPtr,"Unscattered Photon in Lab frame: %e, %e, %e,%e\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3), (ph->r0), (ph->r1), (ph->r2), *(s+0), *(s+1), *(s+2), *(s+3));
            fflush(fPtr);
            fprintf(fPtr,"Fluid Beta: %e, %e, %e\n", *(fluid_beta+0),*(fluid_beta+1), *(fluid_beta+2));
            fflush(fPtr);
            }
            
            fprintf(fPtr,"Old: %e, %e, %e,%e\n", ph->p0, ph->p1, ph->p2, ph->p3);
            fflush(fPtr);
             
            if (((ph+ph_index)->type) == COMPTONIZED_PHOTON)
            {
                fprintf(fPtr, "Before Scattering, In Comov_frame:\n");
                fflush(fPtr);
                fprintf(fPtr, "ph_comov: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                fflush(fPtr);
            }
             */

        
            //then rotate the stokes plane by some angle such that we are in the stokes coordinat eystsem after the lorentz boost
            #if STOKES_SWITCH == ON
            {

                stokesRotation(fluid_beta, (ph_p+1), (ph_p_comov+1), s, fPtr);
                
            }
            #endif
            
            //exit(0);
            //second we generate a thermal electron at the correct temperature
            singleElectron(el_p_comov, fluid_temp, ph_p_comov, rand, fPtr);
            
            /*
            if (((ph+ph_index)->type) == COMPTONIZED_PHOTON)
            {
                fprintf(fPtr,"el_comov: %e, %e, %e,%e\n", *(el_p_comov+0), *(el_p_comov+1), *(el_p_comov+2), *(el_p_comov+3));
                fflush(fPtr);
            }
             */
    
            //third we perform the scattering and save scattered photon 4 monetum in ph_p_comov @ end of function
            event_did_occur=singleScatter(el_p_comov, ph_p_comov, s, rand, fPtr);
        
            /*
            if (((ph+ph_index)->type) == COMPTONIZED_PHOTON)
            {
                fprintf(fPtr,"After Scattering, After Lorentz Boost to Comov frame: %e, %e, %e,%e\n", *(ph_p_comov+0), *(ph_p_comov+1), *(ph_p_comov+2), *(ph_p_comov+3));
                fflush(fPtr);
            }
            */
            
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
                
                /*
                if (((ph+ph_index)->type) == COMPTONIZED_PHOTON)
                {
                    fprintf(fPtr,"Scattered Photon in Lab frame: %e, %e, %e,%e\n\n", *(ph_p+0), *(ph_p+1), *(ph_p+2), *(ph_p+3));
                    fflush(fPtr);
                }
                 */
                
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
                    //energy greater than 1e4 keV
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
        factor=K_B*temp/(M_EL*C_LIGHT*C_LIGHT);
        y_dum=1; //initalize loop to get a random gamma from the distribution of electron velocities
        f_x_dum=0;
        while ((isnan(f_x_dum) !=0) || (y_dum>f_x_dum) )
        {
            
            x_dum=gsl_rng_uniform_pos(rand)*(1+100*factor);
            beta_x_dum=sqrt(1-(1/(x_dum*x_dum)));
            y_dum=gsl_rng_uniform(rand)/2.0;
            
            f_x_dum=x_dum*x_dum*(beta_x_dum/gsl_sf_bessel_Kn (2, 1.0/factor))*exp(-1*x_dum/factor); //
            //fprintf(fPtr,"Choosing a Gamma: xdum: %e, f_x_dum: %e, y_dum: %e\n", x_dum, f_x_dum, y_dum);
        }
        gamma=x_dum;
        
    }
    else
    {

        //printf("In else\n");
        factor=sqrt(K_B*temp/M_EL);
        //calculate a random gamma from 3 random velocities drawn from a gaussian distribution with std deviation of "factor"
        gamma=1.0/sqrt( 1- (pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+ pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)  )); //each vel direction is normal distribution -> maxwellian when multiplied
    }
    
    //fprintf(fPtr,"Chosen Gamma: %e\n",gamma);
    
    beta=sqrt( 1- (1/(gamma*gamma)) );
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
    ph_theta=atan2(sqrt( pow(*(ph_p+2),2)+  pow(*(ph_p+3),2)) , (*(ph_p+1)) );
    
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

double averagePhotonEnergy(struct photon *ph, int num_ph)
{
    //to calculate weighted photon energy in ergs
    int i=0;
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    double e_sum=0, w_sum=0;
    
    #pragma omp parallel for reduction(+:e_sum) reduction(+:w_sum)
    for (i=0;i<num_ph;i++)
    {
        #if CYCLOSYNCHROTRON_SWITCH == ON
        if (((ph+i)->weight != 0)) //dont want account for null or absorbed UNABSORBED_CS_PHOTON photons
        #endif
        {
            e_sum+=(((ph+i)->p0)*((ph+i)->weight));
            w_sum+=((ph+i)->weight);
        }
    }
    
    return (e_sum*C_LIGHT)/w_sum;
}

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg, double *r_avg, FILE *fPtr  )
{
    int temp_max=0, temp_min=INT_MAX,  i=0, count=0, count_synch=0, count_comp=0, count_i=0;
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    double sum=0, avg_r_sum=0, avg_r_sum_synch=0, avg_r_sum_comp=0, avg_r_sum_inject=0;
    
    //printf("Num threads: %d", num_thread);
#pragma omp parallel for num_threads(num_thread) reduction(min:temp_min) reduction(max:temp_max) reduction(+:sum) reduction(+:avg_r_sum) reduction(+:count)
    for (i=0;i<ph_num;i++)
    {
        #if CYCLOSYNCHROTRON_SWITCH == ON
        if (((ph+i)->weight != 0)) //dont want account for null or absorbed UNABSORBED_CS_PHOTON photons
        #endif
        {
            sum+=((ph+i)->num_scatt);
            avg_r_sum+=sqrt(((ph+i)->r0)*((ph+i)->r0) + ((ph+i)->r1)*((ph+i)->r1) + ((ph+i)->r2)*((ph+i)->r2));
            
            //printf("%d %c  %e %e %e %e %e %e\n", i, (ph+i)->type, (ph+i)->p0, (ph+i)->comv_p0, (ph+i)->r0, (ph+i)->r1, (ph+i)->r2, (ph+i)->num_scatt);
            
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
            
            if (((ph+i)->type) == INJECTED_PHOTON )
            {
                avg_r_sum_inject+=sqrt(((ph+i)->r0)*((ph+i)->r0) + ((ph+i)->r1)*((ph+i)->r1) + ((ph+i)->r2)*((ph+i)->r2));
                count_i++;
            }
            
            #if CYCLOSYNCHROTRON_SWITCH == ON
            if ((((ph+i)->type) == COMPTONIZED_PHOTON) || (((ph+i)->type) == UNABSORBED_CS_PHOTON))
            {
                avg_r_sum_comp+=sqrt(((ph+i)->r0)*((ph+i)->r0) + ((ph+i)->r1)*((ph+i)->r1) + ((ph+i)->r2)*((ph+i)->r2));
                count_comp++;
            }
            #endif
            
            count++;
        }
        
        #if CYCLOSYNCHROTRON_SWITCH == ON
        if (((ph+i)->type) == CS_POOL_PHOTON )
        {
            avg_r_sum_synch+=sqrt(((ph+i)->r0)*((ph+i)->r0) + ((ph+i)->r1)*((ph+i)->r1) + ((ph+i)->r2)*((ph+i)->r2));
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

void cylindricalPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    double  gamma_infinity=100, t_comov=1e5, ddensity=3e-7;// the comoving temperature in Kelvin, and the comoving density in g/cm^2
    int i=0;
    double vel=sqrt(1-pow(gamma_infinity, -2.0)), lab_dens=gamma_infinity*ddensity;
    
    fprintf(fPtr, "The Cylindrical Outflow values are: Gamma_infinity=%e, T_comv=%e K, comv dens=%e g/cm^3 \n", gamma_infinity, t_comov, ddensity);
    fflush(fPtr);
    
    for (i=0; i<hydro_data->num_elements; i++)
    {
        ((hydro_data->gamma))[i]=gamma_infinity;
        ((hydro_data->dens))[i]=ddensity;
        ((hydro_data->dens_lab))[i]=lab_dens;
        ((hydro_data->pres))[i]=(A_RAD*pow(t_comov, 4.0))/(3);
        ((hydro_data->temp))[i]=t_comov; //just assign t_comov
        
        
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
            
            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                ((hydro_data->v0))[i]=0;
                ((hydro_data->v1))[i]=vel; //geometry dependent want this to be parallel to jet axis
            #endif

            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel*cos(((hydro_data->r1))[i]);//rhat
                ((hydro_data->v1))[i]=-vel*sin(((hydro_data->r1))[i]);//theta hat direction
            #endif
        
            #if DIMENSIONS == TWO_POINT_FIVE
                //have to make sure that the 3rd vctro direction is set to 0 in 2.5D case
                ((hydro_data->v2))[i]=0;
            #endif
            
        #else

            #if GEOMETRY == CARTESIAN
                ((hydro_data->v0))[i]=0;
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=vel;
            #endif


            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel*cos(((hydro_data->r1))[i]);//rhat
                ((hydro_data->v1))[i]=-vel*sin(((hydro_data->r1))[i]);//theta hat direction
                ((hydro_data->v2))[i]=0;
            #endif

            #if GEOMETRY == POLAR
                ((hydro_data->v0))[i]=0;
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=vel;
            #endif

        #endif
        

    }
    
}

void sphericalPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    double  gamma_infinity=100, lumi=1e54, r00=1e8; //shopuld be 10^57
    //double  gamma_infinity=5, lumi=1e52, r00=1e8; //shopuld be 10^57
    double vel=0, r=0;
    int i=0;
    
    fprintf(fPtr, "The Spherical Outflow values are: Gamma_infinity=%e, Luminosity=%e erg/s, r_0=%e cm \n", gamma_infinity, lumi, r00);
    fflush(fPtr);
    
    for (i=0; i<hydro_data->num_elements; i++)
    {
        if (((hydro_data->r))[i] >= (r00*gamma_infinity))
        {
            ((hydro_data->gamma))[i]=gamma_infinity;
            ((hydro_data->pres))[i]=(lumi*pow(r00, 2.0/3.0)*pow(((hydro_data->r))[i], -8.0/3.0) )/(12.0*M_PI*C_LIGHT*pow(gamma_infinity, 4.0/3.0));
        }
        else
        {
            ((hydro_data->gamma))[i]=((hydro_data->r))[i]/r00;
            ((hydro_data->pres))[i]=(lumi*pow(r00, 2.0))/(12.0*M_PI*C_LIGHT*pow(((hydro_data->r))[i], 4.0) );
        }
        
        ((hydro_data->dens))[i]=lumi/(4*M_PI*pow(((hydro_data->r))[i], 2.0)*pow(C_LIGHT, 3.0)*gamma_infinity*(((hydro_data->gamma))[i]));
        ((hydro_data->dens_lab))[i]=(((hydro_data->dens))[i])*(((hydro_data->gamma))[i]);
        ((hydro_data->temp))[i]=pow(3*(((hydro_data->pres))[i])/(A_RAD) ,1.0/4.0);
        
        vel=sqrt(1-(pow(((hydro_data->gamma))[i], -2.0)));

        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
            
            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=(vel*(((hydro_data->r1))[i]))/r; //geometry dependent want this to be radial
            #endif

            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel;//rhat
                ((hydro_data->v1))[i]=0;//theta hat direction
            #endif
        
            #if DIMENSIONS == TWO_POINT_FIVE
                //have to make sure that the 3rd vctro direction is set to 0 in 2.5D case
                ((hydro_data->v2))[i]=0;
            #endif
            
        #else

            #if GEOMETRY == CARTESIAN
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2)+pow(((hydro_data->r2))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=(vel*(((hydro_data->r1))[i]))/r; //geometry dependent want this to be radial
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif


            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel;//rhat
                ((hydro_data->v1))[i]=0;//theta hat direction
                ((hydro_data->v2))[i]=0;
            #endif

            #if GEOMETRY == POLAR
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r2))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r; //need to figure this out
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif

        #endif

        //fprintf(fPtr,"Gamma: %lf\nR: %lf\nPres: %e\nvel %lf\nX: %lf\nY %lf\nVx: %lf\nVy: %lf\nDens: %e\nLab_Dens: %e\nTemp: %lf\n", *(gamma+i), *(r+i), *(pres+i), vel, *(x+i), *(y+i), *(vx+i), *(vy+i), *(dens+i), *(dens_lab+i), *(temp+i));
    }
    
}

void structuredFireballPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    //This model is provided by Lundman, Peer, Ryde 2014, use this to compare our MCRaT polarization to their polarizations
    double  gamma_0=100, lumi=1e52, r00=1e8, theta_j=1e-2, p=4; //theta_j in paper is 1e-2, 3e-2, 1e-1 and p is 1,2,4
    double T_0=pow(lumi/(4*M_PI*r00*r00*A_RAD*C_LIGHT), 1.0/4.0);
    double eta=0, r_sat=0, r;
    double vel=0, theta_ratio=0;
    int i=0;
    
    fprintf(fPtr, "The Structured Spherical Outflow values are: Gamma_0=%e, Luminosity=%e erg/s, r_0=%e cm, theta_j=%e rad, p=%e \n", gamma_0, lumi, r00, theta_j, p);
    fflush(fPtr);
    
    for (i=0; i<hydro_data->num_elements; i++)
    {
        
        
        theta_ratio=((hydro_data->theta)[i])/theta_j;
        eta=gamma_0/sqrt(1+pow(theta_ratio, 2*p));
        
        if ((hydro_data->theta)[i] >= theta_j*pow(gamma_0/2, 1.0/p))
        {
            //*(gamma+i)=2; //outside with of shear layer have gamma be 2 like in paper
            eta=2.0;
        }
        
        r_sat=eta*r00;
        
        if (((hydro_data->r)[i]) >= r_sat)
        {
            (hydro_data->gamma)[i]=eta;
            (hydro_data->temp)[i]=T_0*pow(r_sat/((hydro_data->r)[i]), 2.0/3.0)/eta;
        }
        else
        {
            (hydro_data->gamma)[i]=((hydro_data->r)[i])/r_sat; //not sure if this is right but it shouldn't matter since we're injecting our photons far from r00
            (hydro_data->temp)[i]=T_0;
        }
        
        vel=sqrt(1-(pow((hydro_data->gamma)[i], -2.0)));
        (hydro_data->dens)[i] = M_P*lumi/(4*M_PI*M_P*C_LIGHT*C_LIGHT*C_LIGHT*eta*vel*((hydro_data->gamma)[i])*((hydro_data->r)[i])*((hydro_data->r)[i])); //equation paper has extra c, but then units dont work out
        (hydro_data->dens_lab)[i]=((hydro_data->dens)[i])*((hydro_data->gamma)[i]);
        (hydro_data->pres)[i]=(A_RAD*pow((hydro_data->temp)[i], 4.0))/(3);
        
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
            
            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=(vel*(((hydro_data->r1))[i]))/r; //geometry dependent want this to be radial
            #endif

            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel;//rhat
                ((hydro_data->v1))[i]=0;//theta hat direction
            #endif
        
            #if DIMENSIONS == TWO_POINT_FIVE
                //have to make sure that the 3rd vctro direction is set to 0 in 2.5D case
                ((hydro_data->v2))[i]=0;
            #endif

        #else

            #if GEOMETRY == CARTESIAN
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2)+pow(((hydro_data->r2))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=(vel*(((hydro_data->r1))[i]))/r; //geometry dependent want this to be radial
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif

            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel;//rhat
                ((hydro_data->v1))[i]=0;//theta hat direction
                ((hydro_data->v2))[i]=0;
            #endif

            #if GEOMETRY == POLAR
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r2))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif

        #endif

        //fprintf(fPtr,"eta: %lf\nr_sat: %lf\nGamma: %lf\nR: %lf\nTheta: %lf\nPres: %e\nvel %lf\nX: %lf\nY %lf\nVx: %lf\nVy: %lf\nDens: %e\nLab_Dens: %e\nTemp: %lf\n\n", eta, r_sat, *(gamma+i), *(r+i), (*(theta+i)), *(pres+i), vel, *(x+i), *(y+i), *(vx+i), *(vy+i), *(dens+i), *(dens_lab+i), *(temp+i));
        
    }
    
}

void phMinMax(struct photon *ph, int ph_num, double *min, double *max, double *min_theta, double *max_theta, FILE *fPtr)
{
    double temp_r_max=0, temp_r_min=DBL_MAX, temp_theta_max=0, temp_theta_min=DBL_MAX;
    int i=0;
    #if defined(_OPENMP)
    int num_thread=omp_get_num_threads();
    #endif
    double ph_r=0, ph_theta=0;
    
#pragma omp parallel for num_threads(num_thread) firstprivate(ph_r, ph_theta) reduction(min:temp_r_min) reduction(max:temp_r_max) reduction(min:temp_theta_min) reduction(max:temp_theta_max)
    for (i=0;i<ph_num;i++)
    {
        if ((ph+i)->weight != 0)
        {
            ph_r=sqrt(((ph+i)->r0)*((ph+i)->r0) + ((ph+i)->r1)*((ph+i)->r1) + ((ph+i)->r2)*((ph+i)->r2));
            ph_theta=acos(((ph+i)->r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            if (ph_r > temp_r_max )
            {
                temp_r_max=ph_r;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            //if ((i==0) || (ph_r<temp_r_min))
            if (ph_r<temp_r_min)
            {
                temp_r_min=ph_r;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
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
        }
    }
    
    *max=temp_r_max;
    *min=temp_r_min;
    *max_theta=temp_theta_max;
    *min_theta=temp_theta_min;
}
