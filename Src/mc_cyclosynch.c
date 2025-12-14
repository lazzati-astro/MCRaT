/*
This file is for the different functions for emitting and absorbing synchrotron photons
*/

#include "mcrat.h"


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
    #if B_FIELD_CALC == TOTAL_E || B_FIELD_CALC == INTERNAL_E
        #if B_FIELD_CALC == INTERNAL_E
            //printf(">> in INTERNAL_E, EPSILON_B is %e \n", EPSILON_B);

            return sqrt(EPSILON_B*8*M_PI*3*el_dens*K_B*temp/2);
        #else B_FIELD_CALC == TOTAL_E
            //printf(">> in TOTAL_E, EPSILON_B is %e \n", EPSILON_B);

            //otherwise calculate B from the total energy
            return sqrt(8*M_PI*EPSILON_B*(el_dens*M_P*C_LIGHT*C_LIGHT+4*A_RAD*temp*temp*temp*temp/3));
        #endif
    #else
        #if B_FIELD_CALC != SIMULATION
            #error B_FIELD_CALC not defined
        #else
            return 0;
        #endif
    #endif
}

double getMagneticFieldMagnitude(struct hydro_dataframe *hydro_data, int hydro_grid_index)
{
    double result=0;
    #if B_FIELD_CALC == TOTAL_E || B_FIELD_CALC == INTERNAL_E
        double el_dens= ((hydro_data->dens)[hydro_grid_index])/M_P;
        result=calcB(el_dens,(hydro_data->temp)[hydro_grid_index]);
    #else
        #if DIMENSIONS == TWO
            result=vectorMagnitude((hydro_data->B0)[hydro_grid_index], (hydro_data->B1)[hydro_grid_index], 0);
        #else
            result=vectorMagnitude((hydro_data->B0)[hydro_grid_index], (hydro_data->B1)[hydro_grid_index], (hydro_data->B2)[hydro_grid_index]);
        #endif
    #endif
    return result;
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

/* Helper: Calculate photon spherical coordinates from Cartesian position */
//TODO: refactor to mcratCoordinateToSphericalCoordinate
static void calculate_photon_position(const struct photon *ph, double *r, double *theta, double *phi)
{
    double x = ph->r0;
    double y = ph->r1;
    double z = ph->r2;
    
    *r = sqrt(x*x + y*y + z*z);
    
    if (*r < DBL_MIN)
    {
        *theta = 0.0;
        *phi = 0.0;
    }
    else
    {
        *theta = acos(z / *r);
        #if DIMENSIONS == THREE
                
                double phi_rad = atan2(y, x);
                *phi = fmod(phi_rad * RAD_TO_DEG + 360.0, 360.0);
        #else
                *phi=0;
        #endif
    }
}

/* Helper: Collect photon statistics and determine ranges */
static int collect_photon_statistics(const struct photonList *photon_list, struct PhotonRangeInfo *info, FILE *fPtr)
{
    *info = (struct PhotonRangeInfo){
        .p0_min = DBL_MAX,
        .p0_max = 0.0,
        .theta_min = DBL_MAX,
        .theta_max = 0.0,
        .phi_min = DBL_MAX,
        .phi_max = 0.0
    };
    
    for (int i = 0; i < photon_list->list_capacity; i++)
    {
        const struct photon *ph = getPhoton(photon_list, i);
                
        
        if ((ph->type != NULL_PHOTON) && (ph->type != CS_POOL_PHOTON))
        {
            
            if (ph->p0 > 0)
            {
                info->p0_min = fmin(info->p0_min, ph->p0);
                info->p0_max = fmax(info->p0_max, ph->p0);
                info->valid_photon_count++;
            }
            
            double r, theta, phi = 0.0;
            calculate_photon_position(ph, &r, &theta, &phi);
            
            info->theta_min = fmin(info->theta_min, theta);
            info->theta_max = fmax(info->theta_max, theta);
            
            #if DIMENSIONS == THREE
            info->phi_min = fmin(info->phi_min, phi);
            info->phi_max = fmax(info->phi_max, phi);
            #endif
        }
        
        if (ph->type == CS_POOL_PHOTON)
        {
            info->synch_photon_count++;
        }
    }
    
    info->log_p0_min = (info->p0_min > 0 && info->p0_max > 0) ? log10(info->p0_min) : 0.0;
    info->log_p0_max = (info->p0_min > 0 && info->p0_max > 0) ? log10(info->p0_max) : 1.0;
    
    return (info->valid_photon_count > 0) ? 1 : 0;
}

/* Helper: Calculate binning parameters based on photon ranges */
static struct BinningParams calculate_binning_params(const struct PhotonRangeInfo *info, int max_photons)
{
    struct BinningParams params = {0};
    params.num_bins = (int)(CYCLOSYNCHROTRON_REBIN_E_PERC * max_photons);
    
    //the size of the bin that we want to produce for spatial binning in theta
    params.num_bins_theta = ceil((info->theta_max - info->theta_min) / (CYCLOSYNCHROTRON_REBIN_ANG * DEG_TO_RAD));
    params.num_bins_phi = 1;
    params.num_avg = 12;
    
    #if DIMENSIONS == THREE
        //the size of the bin that we want to produce for spatial binning in phi, this one is in degrees
        params.num_bins_phi = ceil((info->phi_max - info->phi_min) / CYCLOSYNCHROTRON_REBIN_ANG_PHI);
        params.num_avg = 13;
    #endif
    
    params.total_bins = params.num_bins_theta * params.num_bins;
    #if DIMENSIONS == THREE
        params.total_bins *= params.num_bins_phi;
    #endif
    
    return params;
}

/* Helper: Allocate and initialize all GSL histograms */
static int allocate_histograms(gsl_histogram2d **h_energy_theta, gsl_histogram2d **h_energy_phi, gsl_histogram2d **h_theta_phi, const struct BinningParams *params, const struct PhotonRangeInfo *info, FILE *fPtr)
{
    if (params->num_bins <= 0 || params->num_bins_theta <= 0 || params->num_bins_phi <= 0) {
        fprintf(fPtr, "ERROR: Invalid histogram dimensions\n");
        return 0;
    }
    
    *h_energy_theta = gsl_histogram2d_alloc(params->num_bins, params->num_bins_theta);
    if (!*h_energy_theta) {
        fprintf(fPtr, "ERROR: Failed to allocate energy-theta histogram\n");
        return 0;
    }
    
    double e_eps = (info->log_p0_max - info->log_p0_min) * 1e-6;
    double t_eps = (info->theta_max - info->theta_min) * 1e-6;
    
    gsl_histogram2d_set_ranges_uniform(*h_energy_theta, info->log_p0_min, info->log_p0_max + e_eps, info->theta_min, info->theta_max + t_eps);
    
    #if DIMENSIONS == THREE
        *h_energy_phi = gsl_histogram2d_alloc(params->num_bins, params->num_bins_phi);
        *h_theta_phi = gsl_histogram2d_alloc(params->num_bins_theta, params->num_bins_phi);
        
        if (!*h_energy_phi || !*h_theta_phi)
        {
            fprintf(fPtr, "ERROR: Failed to allocate 3D histograms\n");
            gsl_histogram2d_free(*h_energy_theta);
            if (*h_energy_phi) gsl_histogram2d_free(*h_energy_phi);
            if (*h_theta_phi) gsl_histogram2d_free(*h_theta_phi);
            *h_energy_theta = NULL;
            *h_energy_phi = NULL;
            *h_theta_phi = NULL;
            return 0;
        }
        
        double p_eps = (info->phi_max - info->phi_min) * 1e-6;
        
        gsl_histogram2d_set_ranges_uniform(*h_energy_phi, info->log_p0_min, info->log_p0_max + e_eps, info->phi_min, info->phi_max + p_eps);
        gsl_histogram2d_set_ranges_uniform(*h_theta_phi, info->theta_min, info->theta_max + t_eps, info->phi_min, info->phi_max + p_eps);
    #else
        *h_energy_phi = NULL;
        *h_theta_phi = NULL;
    #endif
    
    return 1;
}

/* Helper: Free all allocated histograms */
static void free_histograms(gsl_histogram2d *h_energy_theta, gsl_histogram2d *h_energy_phi, gsl_histogram2d *h_theta_phi)
{
    if (h_energy_theta)
    {
        gsl_histogram2d_free(h_energy_theta);
    }
    #if DIMENSIONS == THREE
        if (h_energy_phi)
        {
            gsl_histogram2d_free(h_energy_phi);
        }
        if (h_theta_phi)
        {
            gsl_histogram2d_free(h_theta_phi);
        }
    #endif
}

/* Helper: Allocate bin statistics array */
static struct BinStats* allocate_bin_stats(int total_bins, FILE *fPtr)
{
    struct BinStats *stats = calloc(total_bins, sizeof(struct BinStats));
    if (!stats)
    {
        fprintf(fPtr, "ERROR: Failed to allocate %d bin stats\n", total_bins);
    }
    return stats;
}

/* Helper: Free bin statistics array */
static void free_bin_stats(struct BinStats *stats)
{
    free(stats);
}

/* Helper: Calculate safe bin index with bounds checking */
static int calculate_bin_index(int x, int y, int z, const struct BinningParams *params)
{
    if (x < 0 || x >= params->num_bins || y < 0 || y >= params->num_bins_theta)
    {
        return -1;
    }
    
    #if DIMENSIONS == THREE
        if (z < 0 || z >= params->num_bins_phi) return -1;
        return z * params->num_bins * params->num_bins_theta + x * params->num_bins_theta + y;
    #else
        return x * params->num_bins_theta + y;
    #endif
}

/* Helper: Accumulate weighted statistics for each bin */
static int accumulate_bin_statistics(const struct photonList *photon_list, struct BinStats *stats, gsl_histogram2d *h_energy_theta, gsl_histogram2d *h_energy_phi, gsl_histogram2d *h_theta_phi, const struct BinningParams *params, FILE *fPtr)
{
    for (int i = 0; i < photon_list->list_capacity; i++)
    {
        const struct photon *ph = getPhoton(photon_list, i);
                
        if ((ph->type != NULL_PHOTON) && (ph->type != CS_POOL_PHOTON))
        {
            double r, theta, phi = 0.0;
            calculate_photon_position(ph, &r, &theta, &phi);
            
            size_t idx_x = 0, idx_y = 0, idx_z = 0;
            gsl_histogram2d_find(h_energy_theta, log10(ph->p0), theta, &idx_x, &idx_y);
            
            #if DIMENSIONS == THREE
                gsl_histogram2d_find(h_energy_phi, log10(ph->p0), phi, &idx_x, &idx_z);
                gsl_histogram2d_find(h_theta_phi, theta, phi, &idx_y, &idx_z);
            #endif
            
            int bin_idx = calculate_bin_index(idx_x, idx_y, idx_z, params);
            if (bin_idx < 0 || bin_idx >= params->total_bins)
            {
                fprintf(fPtr, "WARNING: Photon %d maps to invalid bin index %d\n", i, bin_idx);
                exit(1);
            }
            
            struct BinStats *s = &stats[bin_idx];
            s->weighted_r += r * ph->weight;
            s->weighted_theta += theta * ph->weight;
            s->weighted_phi_offset += (atan2(ph->p2, ph->p1) - atan2(ph->r1, ph->r0)) * RAD_TO_DEG * ph->weight;
            
            s->weighted_stokes[0] += ph->s0 * ph->weight;
            s->weighted_stokes[1] += ph->s1 * ph->weight;
            s->weighted_stokes[2] += ph->s2 * ph->weight;
            s->weighted_stokes[3] += ph->s3 * ph->weight;
            
            s->weighted_scatt_count += ph->num_scatt * ph->weight;
            s->total_weight += ph->weight;
            
            double phi_dir = fmod(atan2(ph->p2, ph->p1) * RAD_TO_DEG + 360.0, 360.0);
            double theta_dir = acos(ph->p3 / ph->p0) * RAD_TO_DEG;
            
            s->weighted_phi_dir += phi_dir * ph->weight;
            s->weighted_theta_dir += theta_dir * ph->weight;
            s->weighted_energy += ph->p0 * ph->weight;
            
            #if DIMENSIONS == THREE
                s->weighted_phi_pos += phi * ph->weight;
            #endif
        }
    }
    
    return 1;
}

/* Helper: Create rebinned photons from accumulated statistics */
static int create_rebinned_photons(struct photonList *photon_list, const structBinStats *stats, const struct BinningParams *params, int synch_photon_count, FILE *fPtr)
{
    struct photon *rebin_ph = calloc(params->total_bins, sizeof(struct photon));
    if (!rebin_ph)
    {
        fprintf(fPtr, "ERROR: Failed to allocate rebinned photon array\n");
        return -1;
    }
    
    int num_null_rebin_ph = 0;
    
    for (int i = 0; i < params->total_bins; i++)
    {
        const BinStats *s = &stats[i];
        struct photon *new_ph = &rebin_ph[i];
        
        if (s->total_weight <= 0)
        {
            new_ph->type = NULL_PHOTON;
            new_ph->weight = 0;
            new_ph->nearest_block_index = -1;
            new_ph->recalc_properties = 0;
            num_null_rebin_ph++;
        }
        else
        {
            new_ph->type = COMPTONIZED_PHOTON;
            new_ph->weight = s->total_weight;
            
            /* Calculate average values from weighted sums */
            double avg_energy = s->weighted_energy / s->total_weight;
            double avg_phi_dir = s->weighted_phi_dir / s->total_weight;
            double avg_theta_dir = s->weighted_theta_dir / s->total_weight;
            double avg_r = s->weighted_r / s->total_weight;
            double avg_theta_pos = s->weighted_theta / s->total_weight;
            
            /* Set photon momentum components */
            new_ph->p0 = avg_energy;
            new_ph->p1 = avg_energy * sin(avg_theta_dir * DEG_TO_RAD) * cos(avg_phi_dir * DEG_TO_RAD);
            new_ph->p2 = avg_energy * sin(avg_theta_dir * DEG_TO_RAD) * sin(avg_phi_dir * DEG_TO_RAD);
            new_ph->p3 = avg_energy * cos(avg_theta_dir * DEG_TO_RAD);
            
            /* Initialize comoving frame momenta to zero */
            new_ph->comv_p0 = 0;
            new_ph->comv_p1 = 0;
            new_ph->comv_p2 = 0;
            new_ph->comv_p3 = 0;
            
            /* Calculate position phi based on dimensionality */
            double pos_phi;
            #if DIMENSIONS == THREE
                double avg_phi_pos = s->weighted_phi_pos / s->total_weight;
                pos_phi = avg_phi_pos * DEG_TO_RAD;
            #else
                double avg_phi_offset = s->weighted_phi_offset / s->total_weight;
                pos_phi = (avg_phi_dir - avg_phi_offset) * DEG_TO_RAD;
            #endif
            
            /* Set photon position components */
            new_ph->r0 = avg_r * sin(avg_theta_pos) * cos(pos_phi);
            new_ph->r1 = avg_r * sin(avg_theta_pos) * sin(pos_phi);
            new_ph->r2 = avg_r * cos(avg_theta_pos);
            
            /* Set Stokes parameters */
            new_ph->s0 = s->weighted_stokes[0] / s->total_weight;
            new_ph->s1 = s->weighted_stokes[1] / s->total_weight;
            new_ph->s2 = s->weighted_stokes[2] / s->total_weight;
            new_ph->s3 = s->weighted_stokes[3] / s->total_weight;
            
            /* Set scattering count and other properties */
            new_ph->num_scatt = (int)(s->weighted_scatt_count / s->total_weight + 0.5);
            
            //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
            new_ph->nearest_block_index = 0;
            
            //set to 1 so we are sure that we calculate tau values later on
            new_ph->recalc_properties = 1;
        }
    }
    
    /* First, nullify existing Comptonized and Unabsorbed CS photons */
    for (int i = 0; i < photon_list->list_capacity; i++)
    {
        struct photon *ph = getPhoton(photon_list, i);
        if (ph && (ph->type == UNABSORBED_CS_PHOTON || ph->type == COMPTONIZED_PHOTON))
        {
            setNullPhoton(photon_list, i);
        }
    }
    
    /* Add rebinned photons to list */
    addToPhotonList(photon_list, rebin_ph, params->total_bins);
    
    /* Verify all photons were added */
    if (photon_list->num_photons < params->total_bins)
    {
        fprintf(fPtr, "ERROR: Only added %d of %d rebinned photons\n",
                photon_list->num_photons, params->total_bins);
        free(rebin_ph);
        return -1;
    }
    
    free(rebin_ph);
    return num_null_rebin_ph;
}


int rebinCyclosynchCompPhotons(struct photonList *photon_list, int *num_cyclosynch_ph_emit, int *scatt_cyclosynch_num_ph, int max_photons, double thread_theta_min, double thread_theta_max , gsl_rng * rand, FILE *fPtr)
{
    int i=0, j=0, k=0, count=0, count_x=0, count_y=0, count_z=0, count_c_ph=0, end_count=(*scatt_cyclosynch_num_ph), idx=0, num_thread=1;
    #if defined(_OPENMP)
    num_thread=omp_get_num_threads();
    #endif
    int synch_comp_photon_count=0, synch_photon_count=0, num_avg=12, num_bins=(CYCLOSYNCHROTRON_REBIN_E_PERC)*max_photons; //some factor of the max number of photons that is specified in the mc.par file, num bins is also in test function
    double dtheta_bin=CYCLOSYNCHROTRON_REBIN_ANG*M_PI/180; //the size of the bin that we want to produce for spatial binning in theta
    int num_bins_theta=(thread_theta_max-thread_theta_min)/dtheta_bin;//try this many bins such that we have 0.5 degree resolution, can also try to do adaptive binning with constant SNR
    int num_bins_phi=1;//calculate this based on where photons are located in phi, can also try to do adaptive binning with constant SNR
    #if DIMENSIONS == THREE
        double dphi_bin=CYCLOSYNCHROTRON_REBIN_ANG_PHI; //the size of the bin that we want to produce for spatial binning in phi, this one is in degrees
        double ph_phi=0, temp_phi_max=0, temp_phi_min=DBL_MAX, min_range_phi=0, max_range_phi=0;
    #endif
    double avg_values[12]={0}; //number of averages that'll be taken is given by num_avg in above line
    double p0_min=DBL_MAX, p0_max=0, log_p0_min=0, log_p0_max=0;//look at p0 of photons not by frequency since its just nu=p0*C_LIGHT/PL_CONST
    double rand1=0, rand2=0, phi=0, theta=0;
    double min_range=0, max_range=0, min_range_theta=0, max_range_theta=0, energy=0;
    double ph_r=0, ph_theta=0, temp_theta_max=0, temp_theta_min=DBL_MAX;
    //int *synch_comp_photon_idx=NULL; make this an array b/c had issue with deallocating this memory for some reason
    int synch_comp_photon_idx[*scatt_cyclosynch_num_ph], total_bins=0;
    //struct photon *rebin_ph=malloc(num_bins* sizeof (struct photon ));
    int num_null_rebin_ph=0, num_in_bin=0;
    struct photon *tmp_ph=NULL;
    double *tmp_double=NULL;
    int *tmp_int=NULL;
    double count_weight=0;
    //double** avg_values_2d =NULL;

    fprintf(fPtr, "In the rebin func; num_threads %d scatt_cyclosynch_num_ph %d, num_ph %d\n", num_thread, (*scatt_cyclosynch_num_ph), photon_list->num_photons);
    fflush(fPtr);
    
    int min_idx=0, max_idx=0;
    count=0;
    for (i=0;i<photon_list->list_capacity;i++)
    {
        tmp_ph=getPhoton(photon_list, i);
        
        if ((tmp_ph->weight != 0) && ((tmp_ph->type == COMPTONIZED_PHOTON) || (tmp_ph->type == UNABSORBED_CS_PHOTON)) && (tmp_ph->p0 > 0))
        {
            //see if the photon's nu is larger than nu_max or smaller than nu_min
            if ((tmp_ph->p0< p0_min))
            {
                //dont include any absorbed UNABSORBED_CS_PHOTON photons that have negative P0 values
                p0_min= tmp_ph->p0;
                min_idx=i;
                //fprintf(fPtr, "new p0 min %e\n", (p0_min) );
            }
            
            if (tmp_ph->p0> p0_max)
            {
                p0_max= tmp_ph->p0;
                max_idx=i;
                //fprintf(fPtr, "new p0 max %e\n", (p0_max) );
            }
            
            //look at min and max theta of photons
            ph_r=sqrt((tmp_ph->r0)*(tmp_ph->r0) + (tmp_ph->r1)*(tmp_ph->r1) + (tmp_ph->r2)*(tmp_ph->r2));
            ph_theta=acos((tmp_ph->r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            
            if (ph_theta > temp_theta_max )
            {
                temp_theta_max=ph_theta;
                //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            if (ph_theta<temp_theta_min)
            {
                temp_theta_min=ph_theta;
                //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
            }
            
            #if DIMENSIONS == THREE
                ph_phi=fmod(atan((tmp_ph->r1)/ (tmp_ph->r0))*180/M_PI + 360.0,360.0);//want phi to be from 0 to 360 deg
                if (ph_phi > temp_phi_max )
                {
                    temp_phi_max=ph_phi;
                    //fprintf(fPtr, "The new max is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_max, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
                }
                
                if (ph_phi<temp_phi_min)
                {
                    temp_phi_min=ph_phi;
                    //fprintf(fPtr, "The new min is: %e from photon %d with x: %e y: %e z: %e\n", temp_r_min, i, ((ph+i)->r0), (ph+i)->r1, (ph+i)->r2);
                }
            #endif
            
            // also save the index of these photons because they wil become null later on
            synch_comp_photon_idx[count]=i;
            //fprintf(fPtr, "Save index %d\n", i );
            count++;
            
            if (tmp_ph->type == COMPTONIZED_PHOTON)
            {
                //keep track of the number of COMPTONIZED_PHOTON photons so we can know if the array needs to be increased in size, also take num_null_ph into account in doing this
                count_c_ph+=1;
            }
        }
        else if ((tmp_ph->type == CS_POOL_PHOTON) && (tmp_ph->weight != 0))
        {
            synch_photon_count++;
        }
    }
    
    
    num_bins_theta=ceil((temp_theta_max-temp_theta_min)/dtheta_bin);
    #if DIMENSIONS == THREE
        num_bins_phi=ceil((temp_phi_max-temp_phi_min)/dphi_bin);
        num_avg=13;
    #endif
    
    fprintf(fPtr, "Rebin: min, max (keV): %e %e log p0 min, max: %e %e idx: %d %d\n", p0_min*C_LIGHT/1.6e-9,p0_max*C_LIGHT/1.6e-9 , log10(p0_min), log10(p0_max), min_idx, max_idx );
    fprintf(fPtr, "Rebin: min, max (theta in deg): %e %e number of bins %d count: %d\n", temp_theta_min*180/M_PI, temp_theta_max*180/M_PI, num_bins_theta, count );
    #if DIMENSIONS == THREE
        fprintf(fPtr, "Rebin: min, max (phi in deg): %e %e number of bins %d count: %d\n", temp_phi_min, temp_phi_max, num_bins_phi, count );
    #endif
    fflush(fPtr);
    
    if (count != end_count)
    {
        end_count=count; //need this for some reason idk why end_count gets off by 1 compared to what it should be
        fprintf(fPtr, "Rebin: not equal to end_count therefore resetting count to be: %d\n", count );
        fflush(fPtr);
    }
    
    #if DIMENSIONS == THREE
        total_bins=num_bins_phi*num_bins_theta*num_bins;
    #else
        total_bins=num_bins_theta*num_bins;
    #endif
    
    if (total_bins>=max_photons)
    {
        fprintf(fPtr, "The number of rebinned photons, %d, is larger than max_photons %d and will not rebin efficiently. Adjust the parameters such that the number of bins in theta and energy are less than the number of photons that will lead to rebinning.\n",  total_bins, max_photons);
        fflush(fPtr);
        
        printf("Rebin: min, max (theta in deg): %e %e number of bins %d count: %d\n", temp_theta_min*180/M_PI, temp_theta_max*180/M_PI, num_bins_theta, count );
        #if DIMENSIONS == THREE
            printf("Rebin: min, max (phi in deg): %e %e number of bins %d count: %d\n", temp_phi_min, temp_phi_max, num_bins_phi, count );
        #endif
        printf( "In angle range: %e-%e: The number of rebinned photons, %d, is larger than max_photons %d and will not rebin efficiently. Adjust the parameters such that the number of bins in theta and energy are less than the number of photons that will lead to rebinning.\n",  thread_theta_min*180/M_PI, thread_theta_max*180/M_PI, total_bins, max_photons);
        exit(1);

    }

        
    struct photon *rebin_ph=malloc(total_bins* sizeof (struct photon ));
    struct photon *synch_ph=malloc(synch_photon_count* sizeof (struct photon ));
    int synch_photon_idx[synch_photon_count];
    /*
    avg_values_2d = (double**)malloc(total_bins * sizeof(double*));
    for (i = 0; i < total_bins; i++)
    {
        avg_values_2d[i] = (double*)malloc((num_avg) * sizeof(double));
    }
        
    for (i = 0; i < total_bins; i++)
    {
        for (count=0;count<num_avg;count++)
        {
            avg_values_2d[i][count]=0;
        }
    }
     */
    double (*avg_values_2d)[num_avg] = calloc(total_bins, sizeof *avg_values_2d);
    
    gsl_histogram2d * h_energy_theta = gsl_histogram2d_alloc (num_bins, num_bins_theta); //x is for energy  and y is for spatial theta, goes from 0 to pi
    gsl_histogram2d_set_ranges_uniform (h_energy_theta, log10(p0_min), log10(p0_max*(1+1e-6)), temp_theta_min, temp_theta_max*(1+1e-6));
    
    #if DIMENSIONS == THREE
        //make 2D histograms for the other combinations of (phi, theta) and (phi, energy)
        gsl_histogram2d * h_energy_phi = gsl_histogram2d_alloc (num_bins, num_bins_phi); //x is for energy  and y is for spatial theta, goes from 0 to pi
        gsl_histogram2d_set_ranges_uniform (h_energy_phi, log10(p0_min), log10(p0_max*(1+1e-6)), temp_phi_min, temp_phi_max*(1+1e-6));
    
        gsl_histogram2d * h_theta_phi = gsl_histogram2d_alloc (num_bins_theta, num_bins_phi); //x is for energy  and y is for spatial theta, goes from 0 to pi
        gsl_histogram2d_set_ranges_uniform (h_theta_phi, temp_theta_min, temp_theta_max*(1+1e-6), temp_phi_min, temp_phi_max*(1+1e-6));
    #endif


    //populate histogram for photons with nu that falss within the proper histogram bin
    //may not need this loop, can just check if the photon nu falls within the bin edges and do averages etc within next loop
    count=0;
    for (i=0;i<photon_list->list_capacity;i++)
    {
        tmp_ph=getPhoton(photon_list, i);
        
        if ((tmp_ph->weight != 0) && ((tmp_ph->type == COMPTONIZED_PHOTON) || (tmp_ph->type == UNABSORBED_CS_PHOTON)) && (tmp_ph->p0 > 0))
        {
            
            ph_r=sqrt((tmp_ph->r0)*(tmp_ph->r0) + (tmp_ph->r1)*(tmp_ph->r1) + (tmp_ph->r2)*(tmp_ph->r2));
            ph_theta=acos((tmp_ph->r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            
            gsl_histogram2d_increment(h_energy_theta, log10(tmp_ph->p0), ph_theta);
            
            #if DIMENSIONS == THREE
                ph_phi=fmod(atan((tmp_ph->r1)/ (tmp_ph->r0))*180/M_PI + 360.0,360.0);//want phi to be from 0 to 360 deg
                gsl_histogram2d_increment(h_energy_phi, log10(tmp_ph->p0), ph_phi);
                gsl_histogram2d_increment(h_theta_phi, ph_theta, ph_phi);
            #endif
            
            count_weight+=tmp_ph->weight;

        }
        /*
        if ((tmp_ph->type == CS_POOL_PHOTON) && (tmp_ph->weight != 0))
        {
            //save the sych photons here because they may get written over later and corrupted
            (synch_ph+count)->p0=tmp_ph->p0;
            (synch_ph+count)->p1=tmp_ph->p1;
            (synch_ph+count)->p2=tmp_ph->p2;
            (synch_ph+count)->p3=tmp_ph->p3;
            (synch_ph+count)->comv_p0=tmp_ph->comv_p0;
            (synch_ph+count)->comv_p1=tmp_ph->comv_p1;
            (synch_ph+count)->comv_p2=tmp_ph->comv_p2;
            (synch_ph+count)->comv_p3=tmp_ph->comv_p3;
            (synch_ph+count)->r0=tmp_ph->r0;
            (synch_ph+count)->r1= tmp_ph->r1;
            (synch_ph+count)->r2=tmp_ph->r2; //y coordinate in flash becomes z coordinate in MCRaT
            (synch_ph+count)->s0=tmp_ph->s0; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
            (synch_ph+count)->s1=tmp_ph->s1;
            (synch_ph+count)->s2=tmp_ph->s2;
            (synch_ph+count)->s3=tmp_ph->s3;
            (synch_ph+count)->num_scatt=tmp_ph->num_scatt;
            (synch_ph+count)->weight=tmp_ph->weight;
            (synch_ph+count)->nearest_block_index=tmp_ph->nearest_block_index; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
            synch_photon_idx[count]=i;
            count++;
        }
         */
         
    }
    
    //fprintf(fPtr, "counted_weight 1 %e\n", count_weight);
    //fflush(fPtr);
    count_weight=0;
    //gsl_histogram2d_fprintf(fPtr, h_energy_theta, "%g", "%g");
    //have new loop where we loop over the photons, calculate their theta, E, and phi. Have 2D array of size total_bins*12 that we use to calculate averages
    //double avg_values_2d[total_bins][12]={0}; //number of averages that'll be taken is given by num_avg
    
    

        
 
    count=0; //used to claculate where in array we should save averages to
    for (i=0;i<photon_list->list_capacity;i++)
    {
        tmp_ph=getPhoton(photon_list, i);
        
        if ((tmp_ph->weight != 0) && ((tmp_ph->type == COMPTONIZED_PHOTON) || (tmp_ph->type == UNABSORBED_CS_PHOTON)) && (tmp_ph->p0 > 0))
        {
            
            ph_r=sqrt((tmp_ph->r0)*(tmp_ph->r0) + (tmp_ph->r1)*(tmp_ph->r1) + (tmp_ph->r2)*(tmp_ph->r2));
            ph_theta=acos((tmp_ph->r2) /ph_r); //this is the photons theta psition in the FLASH grid, gives in radians
            
            gsl_histogram2d_find(h_energy_theta, log10(tmp_ph->p0), ph_theta, &count_x, &count_y);
            
            #if DIMENSIONS == THREE
                ph_phi=fmod(atan((tmp_ph->r1)/ (tmp_ph->r0))*180/M_PI + 360.0,360.0);//want phi to be from 0 to 360 deg
                gsl_histogram2d_find(h_energy_phi, log10(tmp_ph->p0), ph_phi, &count_x, &count_z);
                gsl_histogram2d_find(h_theta_phi, ph_theta, ph_phi, &count_y, &count_z);
            #endif
            
            //calculate index in avg_values_2d where the histogram bins lie
            count=count_z*num_bins*num_bins_theta+count_x*num_bins_theta+count_y; //count_z=0 if this is 2D case otherwise
            //fprintf(fPtr, "Photon %d with count %d weight %e  and avg_value_weight %e\n", i, count,  tmp_ph->weight, avg_values_2d[count][8]);
            //fflush(fPtr);
            
            //calculate cumulative averages
            avg_values_2d[count][0] += ph_r*tmp_ph->weight; // doing r, theta averages in space
            avg_values_2d[count][1] += ph_theta*tmp_ph->weight;
            avg_values_2d[count][2] += ((atan(tmp_ph->p2/(tmp_ph->p1))*180/M_PI)-(atan((tmp_ph->r1)/ (tmp_ph->r0))*180/M_PI))*tmp_ph->weight;// look at delta \phi between the 4 mometum and its location
            avg_values_2d[count][3] += tmp_ph->s0*tmp_ph->weight;
            avg_values_2d[count][4] += tmp_ph->s1*tmp_ph->weight;
            avg_values_2d[count][5] += tmp_ph->s2*tmp_ph->weight;
            avg_values_2d[count][6] += tmp_ph->s3*tmp_ph->weight;
            avg_values_2d[count][7] += tmp_ph->num_scatt*tmp_ph->weight;
            avg_values_2d[count][8] += tmp_ph->weight;
            
            //get avg theta and phi
            {
                avg_values_2d[count][9] += fmod(atan2(tmp_ph->p2,(tmp_ph->p1))*180/M_PI + 360.0,360.0) *tmp_ph->weight;
                avg_values_2d[count][10] += (180/M_PI)*acos((tmp_ph->p3)/(tmp_ph->p0))*tmp_ph->weight;

            }
            
            avg_values_2d[count][11] +=tmp_ph->p0*tmp_ph->weight;
            #if DIMENSIONS == THREE
                avg_values_2d[count][12] += ph_phi*tmp_ph->weight;
            #endif

        }
    }
    
    //save the avg values to the rebin_arrays
    num_null_rebin_ph=0;
    for (i=0;i<total_bins;i++)
    {
        //test if weight== 0 in 2D array, means there are no photons in that theta,E, phi bin
        if (avg_values_2d[i][8]==0)
        {
            (rebin_ph+i)->type = NULL_PHOTON;
            //fprintf(fPtr, "bin %e-%e has %e photons\n", pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, 0.0);

            
            (rebin_ph+i)->p0=1;
            (rebin_ph+i)->p1=0;
            (rebin_ph+i)->p2=0;
            (rebin_ph+i)->p3=0;
            (rebin_ph+i)->comv_p0=0;
            (rebin_ph+i)->comv_p1=0;
            (rebin_ph+i)->comv_p2=0;
            (rebin_ph+i)->comv_p3=0;
            (rebin_ph+i)->r0=0;
            (rebin_ph+i)->r1= 0;
            (rebin_ph+i)->r2=0;
            (rebin_ph+i)->s0=1; // stokes parameterized are normalized such that I always =1
            (rebin_ph+i)->s1=0;
            (rebin_ph+i)->s2=0;
            (rebin_ph+i)->s3=0;
            (rebin_ph+i)->num_scatt=0;
            (rebin_ph+i)->weight=0;
            (rebin_ph+i)->nearest_block_index=-1; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
            (rebin_ph+i)->recalc_properties=0; //set to 1 so we are sure that we calculate tau values later on

            count_weight+=(rebin_ph+i)->weight;
            num_null_rebin_ph++;

        }
        else
        {
            energy=avg_values_2d[i][11]/avg_values_2d[i][8];
            
            phi=avg_values_2d[i][9]/avg_values_2d[i][8];
            theta=avg_values_2d[i][10]/avg_values_2d[i][8];
            
            (rebin_ph+i)->type = COMPTONIZED_PHOTON;
                            
            
            (rebin_ph+i)->p0=energy;
            (rebin_ph+i)->p1=energy*sin(theta*M_PI/180)*cos(phi*M_PI/180);
            (rebin_ph+i)->p2=energy*sin(theta*M_PI/180)*sin(phi*M_PI/180);
            (rebin_ph+i)->p3=energy*cos(theta*M_PI/180);
            (rebin_ph+i)->comv_p0=0;
            (rebin_ph+i)->comv_p1=0;
            (rebin_ph+i)->comv_p2=0;
            (rebin_ph+i)->comv_p3=0;
            
            #if DIMENSIONS == THREE
                //actually use average phi here since these photons are already in teh same phi bin
                rand1=(M_PI/180)*(avg_values_2d[i][12]/avg_values_2d[i][8]);
            #else
                //calculate the rebinned photon's positional phi as a displacement from its 4-momentum phi direction in 2D, dont care where photons are in phi because of cylilndrical symmetry
                rand1=(M_PI/180)*(phi-avg_values_2d[i][2]/avg_values_2d[i][8]);
            #endif

            
            
            (rebin_ph+i)->r0= (avg_values_2d[i][0]/avg_values_2d[i][8])*sin(avg_values_2d[i][1]/avg_values_2d[i][8])*cos(rand1); //avg_values[0]/avg_values[8]; now do avg r * avg theta * random phi
            (rebin_ph+i)->r1= (avg_values_2d[i][0]/avg_values_2d[i][8])*sin(avg_values_2d[i][1]/avg_values_2d[i][8])*sin(rand1); //avg_values[1]/avg_values[8];
            (rebin_ph+i)->r2= (avg_values_2d[i][0]/avg_values_2d[i][8])*cos(avg_values_2d[i][1]/avg_values_2d[i][8]); //avg_values[2]/avg_values[8];

            (rebin_ph+i)->s0=avg_values_2d[i][3]/avg_values_2d[i][8]; // stokes parameterized are normalized such that I always =1
            (rebin_ph+i)->s1=avg_values_2d[i][4]/avg_values_2d[i][8];
            (rebin_ph+i)->s2=avg_values_2d[i][5]/avg_values_2d[i][8];
            (rebin_ph+i)->s3=avg_values_2d[i][6]/avg_values_2d[i][8];
            (rebin_ph+i)->num_scatt=avg_values_2d[i][7]/avg_values_2d[i][8];
            (rebin_ph+i)->weight=avg_values_2d[i][8];
            (rebin_ph+i)->nearest_block_index=0; //hopefully this is not actually the block that this photon's located in b/c we need to get the 4 mometum in the findNearestProperties function
            (rebin_ph+i)->recalc_properties=1; //set to 1 so we are sure that we calculate tau values later on

            //fprintf(fPtr, "bin %e-%e, %e-%e has %e photons: Theta of averages photon is: %e\n",pow(10, min_range)*C_LIGHT/1.6e-9, pow(10,max_range)*C_LIGHT/1.6e-9, min_range_theta, max_range_theta, gsl_histogram2d_get(h_energy_theta, count_x, count_y), ph_theta);
            //fflush(fPtr);
            
            count_weight+=(rebin_ph+i)->weight;

        }
    }
    
    //free the dynamically allocated memory
    /*
    for (i = 0; i < total_bins; i++)
        free(avg_values_2d[i]);
     */
    free(avg_values_2d);
    

    //end second loop that can replace original nested loop

    
    
    //fprintf(fPtr, "counted_weight 2 %e\n", count_weight);
    //fflush(fPtr);
    
    //exit(0);
    
    //go through all the phtons and ID where the comptonized/unabsorbed CS ones are and set them to be NULL_PHOTONS. These will open up these slots to be written into later. It will also update our counter of num_null_photons in the struct
    for (i=0;i<photon_list->list_capacity;i++)
    {
        tmp_ph=getPhoton(photon_list, i);

        if ((tmp_ph->type == UNABSORBED_CS_PHOTON) || (tmp_ph->type == COMPTONIZED_PHOTON))
        {
            setNullPhoton(photon_list, i);
        }
    }

        
    //add all the CS pool photons (again out of abundance of caution), only if we set the indexes of these photons to be NULL above
    //addToPhotonList(photon_list, synch_ph, synch_photon_count);
        
    //add all the rebinned photons to the list
    addToPhotonList(photon_list, rebin_ph, total_bins);
    
    //make sure that all the rebinned photons have been saved
    /*
    if (count<total_bins)
    {
        fprintf(fPtr, "There was an issue where MCRaT was not able to save all of the rebinned photons\n");
        printf("TThere was an issue where MCRaT was not able to save all of the rebinned photons\n");
        fflush(fPtr);
        exit(1);
    }
     */
    
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

    
    
    *scatt_cyclosynch_num_ph=total_bins-num_null_rebin_ph;
    *num_cyclosynch_ph_emit=total_bins+synch_photon_count-num_null_rebin_ph; //include the emitted synch photons and exclude any of those that are null
    //*num_null_ph=j; //was using j before but i have no idea why its not counting correctly
        
    //fprintf(fPtr, "orig null_ph: %d Calc num_ph: %d counted null_ph: %d  num_null_rebin_ph: %d old scatt_cyclosynch_num_ph: %d new scatt_cyclosynch_num_ph: %d\n", *num_null_ph, (*num_ph), j, num_null_rebin_ph, *scatt_cyclosynch_num_ph, total_bins-num_null_rebin_ph  );
    //fflush(fPtr);
    
    //exit(2);

    
    gsl_histogram2d_free (h_energy_theta);
    #if DIMENSIONS == THREE
        gsl_histogram2d_free(h_energy_phi);
    #endif
    free(rebin_ph);
    free(synch_ph);
    
    return num_null_rebin_ph;
}



int photonEmitCyclosynch(struct photonList *photon_list, double r_inj, double ph_weight, int maximum_photons, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng *rand, int inject_single_switch, int scatt_ph_index, FILE *fPtr)
{
    double rmin=0, rmax=0, max_photons=CYCLOSYNCHROTRON_REBIN_E_PERC*maximum_photons; //have 10% as default, can change later need to figure out how many photons across simulations I want emitted
    double ph_weight_adjusted=0, position_phi=0;
    double dimlesstheta=0, nu_c=0, el_dens=0, error=0, ph_dens_calc=0, max_jnu=0, b_field=0;
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
            #if DIMENSIONS == THREE
                //want inner corner to be close to origin, therfore ned to have abs for 3D cartesian with negative coordinates, shouldnt affect the other geometry systems since theyre all defined from r=0, theta=0, phi=0
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, fabs((hydro_data->r0)[i])-0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])-0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])-0.5*(hydro_data->r2_size)[i]);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, fabs((hydro_data->r0)[i])+0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])+0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])+0.5*(hydro_data->r2_size)[i]);
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
                #if DIMENSIONS == THREE
                    //want inner corner to be close to origin, therfore ned to have abs for 3D cartesian with negative coordinates, shouldnt affect the other geometry systems since theyre all defined from r=0, theta=0, phi=0
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, fabs((hydro_data->r0)[i])-0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])-0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])-0.5*(hydro_data->r2_size)[i]);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, fabs((hydro_data->r0)[i])+0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])+0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])+0.5*(hydro_data->r2_size)[i]);
                #else
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
                #endif

                if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  < rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner < theta_max))
                {
                    //set parameters fro integration fo phtoons spectrum
                    el_dens= ((hydro_data->dens)[i])/M_P;

                    b_field=getMagneticFieldMagnitude(hydro_data,i);
                    nu_c=calcCyclotronFreq(b_field);
                    dimlesstheta=calcDimlessTheta( (hydro_data->temp)[i]);
                    //fprintf(fPtr, "B field is:%e %e %e with magnitude %e nu_c is %e at r=%e\n", (hydro_data->B0)[i], (hydro_data->B1)[i], (hydro_data->B2)[i], b_field, nu_c, (hydro_data->r)[i]);
                    //fflush(fPtr);

                    //printf("Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e\n",*(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta);

                    params[0] = (hydro_data->temp)[i]; //nu_c;
                    params[1] = dimlesstheta;
                    params[2] = el_dens;
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
            else
            {
            
            fprintf(fPtr, "dens: %d, photons: %d, adjusted weight: %e\n", *(ph_dens+(j-1)), ph_tot, ph_weight_adjusted);
            fflush(fPtr);
            }
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
    
    /*FIND OUT WHICH PHOTONS IN ARRAY ARE OLD/WERE ABSORBED AND IDENTIFY THIER INDEXES AND HOW MANY, dont subtract this from ph_tot @ the end, WILL NEED FOR PRINT PHOTONS
    #pragma omp parallel for num_threads(num_thread) reduction(+:null_ph_count)
    for (i=0;i<*num_ph;i++)
    {
        if (((*ph_orig)[i].weight == 0)) //if photons are null COMPTONIZED_PHOTON photons and not absorbed UNABSORBED_CS_PHOTON photons
        {
            null_ph_count+=1;
        }
    }
    */ //not needed since using the photonList struct
    
    //allocate memory for that many photons and also allocate memory to hold comoving 4 momentum of each photon and the velocity of the fluid
    ph_emit=malloc (ph_tot * sizeof (struct photon ));
    p_comv=malloc(4*sizeof(double));
    boost=malloc(4*sizeof(double));
    l_boost=malloc(4*sizeof(double));
    
    
//    if (null_ph_count < ph_tot)
//    {
//        //if the totoal number of photons to be emitted is larger than the number of null phtons curently in the array, then have to grow the array
//        //need to realloc memory to hold the old photon info and the new emitted photon's info
//        //fprintf(fPtr, "Emit: Allocating %d space\n", ((*num_ph)+ph_tot-null_ph_count));
//        //fflush(fPtr);
//        /*
//        tmp=realloc(*ph_orig, ((*num_ph)+ph_tot-null_ph_count)* sizeof (struct photon )); //may have to look into directly doubling (or *1.5) number of photons each time we need to allocate more memory, can do after looking at profiling for "just enough" memory method
//        if (tmp != NULL)
//        {
//            // everything ok
//            *ph_orig = tmp;
//
//        }
//        else
//        {
//            // problems!!!!
//            printf("Error with reserving space to hold old and new photons\n");
//            exit(0);
//        }
//        */
//        //also expand memory of other arrays
//        /* this isnt needed since the time steps are contained in the photon struct
//        tmp_double=realloc(*all_time_steps, ((*num_ph)+ph_tot-null_ph_count)*sizeof(double));
//        if (tmp_double!=NULL)
//        {
//            *all_time_steps=tmp_double;
//        }
//        else
//        {
//            printf("Error with reallocating space to hold data about each photon's time step until an interaction occurs\n");
//        }
//        */
//        tmp_int=realloc(*sorted_indexes, ((*num_ph)+ph_tot-null_ph_count)*sizeof(int));
//        if (tmp_int!=NULL)
//        {
//            *sorted_indexes=tmp_int;
//        }
//        else
//        {
//            printf("Error with reallocating space to hold data about the order in which each photon would have an interaction\n");
//        }
//        
//        net_ph=(ph_tot-null_ph_count);
//        null_ph_count=ph_tot; // use this to set the photons recently allocated as null phtoons (this can help if we decide to directly double (or *1.5) number of photons each time we need to allocate more memory, then use factor*((*num_ph)+ph_tot)-(*num_ph)
//        null_ph_indexes=malloc((ph_tot+null_ph_count)*sizeof(int));
//        j=0;
//        for (i=((*num_ph)+net_ph)-1;i >=0 ;i--)
//        {
//            //fprintf(fPtr, "idx %d\n", i);
//            //fflush(fPtr);
//            if (((*ph_orig)[i].weight == 0)   || (i >= *num_ph))
//            {
//                //preset values for the the newly created spots to hold the emitted phtoons in
//                (*ph_orig)[i].weight=0;
//                (*ph_orig)[i].nearest_block_index=-1;
//                *(null_ph_indexes+j)=i; //save this information so we can use the same syntax for both cases in saving the emitted photon data
//                //fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
//                //fflush(fPtr);
//                j++;
//            }
//        }
//        count_null_indexes=ph_tot; //use this to count the number fo null photons we have actually created, (this can help if we decide to directly double (or *1.5) number of photons each time we need to allocate more memory, then use factor*((*num_ph)+ph_tot)-(*num_ph)
//        
//        //loop through the original set of photons to see if
//        
//        //fprintf(fPtr,"Val %d\n", (*(null_ph_indexes+count_null_indexes-1)));
//        *num_ph+=net_ph; //update number of photons
//        *num_null_ph=ph_tot-null_ph_count; //((*num_ph)+ph_tot)-(*num_ph)-ph_tot; //reserved space - emitted photons-original photons
//        //fprintf(fPtr,"New Num PH %d\nNew null hum_ph %d\n", *num_ph, *num_null_ph);
//        //fflush(fPtr);
//    }
//    else
//    {
//        //otherwise need to find the indexes of these null photons to save the newly emitted photons in them, start searching from the end of the array to efficiently find them
//        //dont need to update the number of photons here
//        null_ph_indexes=malloc(null_ph_count*sizeof(int));
//        j=0;
//        for (i=(*num_ph)-1;i>=0;i--)
//        {
//            if ((*ph_orig)[i].weight == 0)  //if photons are null COMPTONIZED_PHOTON photons and not absorbed UNABSORBED_CS_PHOTON photons
//            {
//                // if the weight is 0, this is a photons that has been absorbed and is now null
//                *(null_ph_indexes+j)=i;
//                j++;
//                //fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
//                //fflush(fPtr);
//                
//                if (j == null_ph_count)
//                {
//                    i=-1; //have found al the indexes and can exit the loop, dont want to do this so we can do the first part of the if statement
//                }
//            }
//            
//        }
//        
//        count_null_indexes=null_ph_count;
//        
//        *num_null_ph=null_ph_count-ph_tot;
//        
//    }
    //the above chunk of code is meant to allocate extra memory or idenify the location of null photons in the photon array. These are no longer needed since using the photonList struct automatically determines if the list needs to be grown and sets null photons to the appropriate values which can be over written below. Now we allocate an array of photons that will be added to the list at the end of this function
    
    //allocate memory for that many photons and also allocate memory to hold comoving 4 momentum of each photon and the velocity of the fluid
    
    if (inject_single_switch == 0)
    {
        //go through blocks and assign random energies/locations to proper number of photons
        net_ph=ph_tot; //save this value to check later on if we have created all the photons that we expect to create
        ph_tot=0;
        for (i=0;i< hydro_data->num_elements;i++)
        {
            #if DIMENSIONS == THREE
                //want inner corner to be close to origin, therfore ned to have abs for 3D cartesian with negative coordinates, shouldnt affect the other geometry systems since theyre all defined from r=0, theta=0, phi=0
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, fabs((hydro_data->r0)[i])-0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])-0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])-0.5*(hydro_data->r2_size)[i]);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, fabs((hydro_data->r0)[i])+0.5*(hydro_data->r0_size)[i], fabs((hydro_data->r1)[i])+0.5*(hydro_data->r1_size)[i], fabs((hydro_data->r2)[i])+0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (hydro_data->r0)[i]-0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]-0.5*(hydro_data->r1_size)[i], 0);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (hydro_data->r0)[i]+0.5*(hydro_data->r0_size)[i], (hydro_data->r1)[i]+0.5*(hydro_data->r1_size)[i], 0);
            #endif

            if ((rmin <= r_grid_outercorner) && (r_grid_innercorner  < rmax ) && (theta_grid_outercorner >= theta_min) && (theta_grid_innercorner < theta_max))
            {
                
                b_field=getMagneticFieldMagnitude(hydro_data,i);
                nu_c=calcCyclotronFreq(b_field);

                dimlesstheta=calcDimlessTheta( (hydro_data->temp)[i]);
                
                for(j=0;j<( *(ph_dens+k) ); j++ )
                {
                    //printf("flash_array_idx: %d Temp %e, el_dens %e, B %e, nu_c %e, dimlesstheta %e\n",i, *(temp+i), el_dens, calcB(el_dens, *(temp+i), epsilon_b), nu_c, dimlesstheta);

                    fr_dum=nu_c; //set the frequency directly to the cyclotron frequency
                    //fprintf(fPtr, "%lf\n ",fr_dum);
                    //exit(0);
                    #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
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
      
                    //printf("%lf, %lf, %lf\n", *(boost+0), *(boost+1), *(boost+2));
                    
                    //boost to lab frame
                    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
                    //printf("Assigning values to struct\n");
                    
                    idx=ph_tot; //(*(null_ph_indexes+count_null_indexes-1));
                    //fprintf(fPtr, "Placing photon in index %d\n", idx);
                    ph_emit[idx].p0=(*(l_boost+0));
                    ph_emit[idx].p1=(*(l_boost+1));
                    ph_emit[idx].p2=(*(l_boost+2));
                    ph_emit[idx].p3=(*(l_boost+3));
                    ph_emit[idx].comv_p0=(*(p_comv+0));
                    ph_emit[idx].comv_p1=(*(p_comv+1));
                    ph_emit[idx].comv_p2=(*(p_comv+2));
                    ph_emit[idx].comv_p3=(*(p_comv+3));
                    
                    #if DIMENSIONS == THREE
                        hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
                    #else
                        hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
                    #endif
                    ph_emit[idx].r0= cartesian_position_rand_array[0]; //put photons @center of the box with random phi
                    ph_emit[idx].r1= cartesian_position_rand_array[1] ;
                    ph_emit[idx].r2= cartesian_position_rand_array[2]; //y coordinate in flash becomes z coordinate in MCRaT
                    
                    //fprintf(fPtr,"%d %e %e %e\n", ph_tot, ph_emit[idx].r0, ph_emit[idx].r1, ph_emit[idx].r2);
                    //fflush(fPtr);
                    
                    ph_emit[idx].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
                    ph_emit[idx].s1=0;
                    ph_emit[idx].s2=0;
                    ph_emit[idx].s3=0;
                    ph_emit[idx].num_scatt=0;
                    ph_emit[idx].weight=ph_weight_adjusted;
                    ph_emit[idx].nearest_block_index=0; //these photons can be scattered
                    ph_emit[idx].type=CS_POOL_PHOTON;
                    ph_emit[idx].recalc_properties=1; //set to 1 so we are sure that we calculate tau values later on

                    //printf("%d\n",ph_tot);
                    ph_tot++; //count how many photons have been emitted
                    //count_null_indexes--; //keep track fo the null photon indexes
                    
                    if (net_ph==ph_tot)
                    {
                        //if we have created all the emitted photons
                        //ph_tot is equal to what it used to be
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
        idx=0; //(*(null_ph_indexes+count_null_indexes-1));
        tmp=getPhoton(photon_list, scatt_ph_index);
        i=tmp->nearest_block_index;

        b_field=getMagneticFieldMagnitude(hydro_data, i);
        nu_c=calcCyclotronFreq(b_field);
        
        fr_dum=nu_c; //_scatt; //set the frequency directly to the cyclotron frequency
        //fprintf(fPtr, "%lf %d\n ",fr_dum, (*ph_orig)[scatt_ph_index].nearest_block_index);
        //exit(0);
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
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
        //printf("%lf, %lf, %lf\n", *(boost+0), *(boost+1), *(boost+2));
        
        //boost to lab frame
        lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);
        
        //fprintf(fPtr, "Placing photon in index %d\n", idx);
        ph_emit[idx].p0=(*(l_boost+0));
        ph_emit[idx].p1=(*(l_boost+1));
        ph_emit[idx].p2=(*(l_boost+2));
        ph_emit[idx].p3=(*(l_boost+3));
        ph_emit[idx].comv_p0=(*(p_comv+0));
        ph_emit[idx].comv_p1=(*(p_comv+1));
        ph_emit[idx].comv_p2=(*(p_comv+2));
        ph_emit[idx].comv_p3=(*(p_comv+3));

        #if DIMENSIONS == THREE
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], (hydro_data->r2)[i]);
        #else
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i], (hydro_data->r1)[i], position_phi);
        #endif
        ph_emit[idx].r0= cartesian_position_rand_array[0]; //put photons @center of the box with random phi
        ph_emit[idx].r1= cartesian_position_rand_array[1] ;
        ph_emit[idx].r2= cartesian_position_rand_array[2]; //y coordinate in flash becomes z coordinate in MCRaT

        ph_emit[idx].s0=1; //initalize stokes parameters as non polarized photon, stokes parameterized are normalized such that I always =1
        ph_emit[idx].s1=0;
        ph_emit[idx].s2=0;
        ph_emit[idx].s3=0;
        ph_emit[idx].num_scatt=0;
        ph_emit[idx].weight=tmp->weight;
        ph_emit[idx].nearest_block_index=i; //these photons can be scattered
        ph_emit[idx].type=CS_POOL_PHOTON;
        ph_emit[idx].recalc_properties=1; //set to 1 so we are sure that we calculate tau values later on

        
        //change position of scattered synchrotron photon to be random in the hydro grid
        position_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r0_size)[i])-((hydro_data->r0_size)[i])/2.0; //choose between -size/2 to size/2
        position2_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r1_size)[i])-((hydro_data->r1_size)[i])/2.0;
        #if DIMENSIONS == THREE
            position3_rand=gsl_rng_uniform_pos(rand)*((hydro_data->r2_size)[i])-((hydro_data->r2_size)[i])/2.0;
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, (hydro_data->r2)[i]+position3_rand);
        #else
            hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array, (hydro_data->r0)[i]+position_rand, (hydro_data->r1)[i]+position2_rand, position_phi);
        #endif

        //assign random position
        tmp->r0=cartesian_position_rand_array[0];
        tmp->r1=cartesian_position_rand_array[1];
        tmp->r2=cartesian_position_rand_array[2];
        
    }
    //printf("At End of function\n");
    //exit(0);
    
    //add the whole array to our photon list struct
    addToPhotonList(photon_list, ph_emit, ph_tot);

    free(null_ph_indexes);
    free(ph_dens); free(p_comv); free(boost); free(l_boost);
    free(ph_emit);
    
    gsl_integration_workspace_free (w);
    
    return ph_tot;
}

double phAbsCyclosynch(struct photonList *photon_list, int *num_abs_ph, int *scatt_cyclosynch_num_ph, struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    int i=0, count=0, abs_ph_count=0, synch_ph_count=0, num_thread=1;
    int other_count=0;
    #if defined(_OPENMP)
    num_thread=omp_get_num_threads();
    #endif

    double el_dens=0, nu_c=0, abs_count=0, b_field=0;
    struct photon *ph=NULL;

    fprintf(fPtr, "In phAbsCyclosynch func begin: abs_ph_count: %d synch_ph_count: %d scatt_cyclosynch_num_ph: %d num_threads: %d\n", abs_ph_count, synch_ph_count, *scatt_cyclosynch_num_ph, num_thread);
    
    *scatt_cyclosynch_num_ph=0;//set thsi equal to 0, to recount in this function and get prepared for the next frame
    
    #pragma omp parallel for num_threads(num_thread) firstprivate(b_field, el_dens, nu_c) reduction(+:abs_ph_count)
    for (i=0;i<photon_list->list_capacity;i++)
    {
        ph=getPhoton(photon_list, i);
        
        if ((ph->weight != 0) && (ph->nearest_block_index != -1))
        {
            // if the photon isnt a null photon already, see if it should be absorbed
            
            b_field=getMagneticFieldMagnitude(hydro_data, ph->nearest_block_index);
            nu_c=calcCyclotronFreq(b_field);

            //printf("photon %d has lab nu %e comv frequency %e and nu_c %e with FLASH grid number %d\n", i, ph->p0*C_LIGHT/PL_CONST, ph->comv_p0*C_LIGHT/PL_CONST, nu_c, ph->nearest_block_index);
            if ((ph->comv_p0*C_LIGHT/PL_CONST <= nu_c) || (ph->type == CS_POOL_PHOTON))
            {
                //if the photon has a frequency less that nu_c, it should be absorbed and becomes a null photon
                //preset values for the the newly created spots to hold the emitted phtoons in;
                
                //if this is a synchrotron photons or photons that have been scattered that were once synch photons in this frame
                //fprintf(fPtr,"photon %d being absorbed\n", i);
                abs_ph_count++;

                if ((ph->type != INJECTED_PHOTON) && (ph->type != UNABSORBED_CS_PHOTON) )
                {
                    if (ph->type == CS_POOL_PHOTON)
                    {
                        synch_ph_count++;
                    }
                }
                else
                {
                    //have an injected photon or UNABSORBED_CS_PHOTON (previous COMPTONIZED_PHOTON photon) that has a nu that can be absorbed
                    abs_count+=ph->weight;
                    ph->p0=-1; //set its energy negative so we know for later analysis that it can't be used and its been absorbed,
                }
                
                setNullPhoton(photon_list, i);

            }
            else
            {
                if ((ph->type == COMPTONIZED_PHOTON) || (ph->type == UNABSORBED_CS_PHOTON) )
                {
                    //if the photon is a COMPTONIZED_PHOTON phton (scattered synch photon from the current frame) or a UNABSORBED_CS_PHOTON photon (scattered synch photon) from an old frame
                    //count how many of these there are
                    *scatt_cyclosynch_num_ph+=1;
                }

            }
            /*
            else
            {
                //if the phootn isnt going to be absorbed, see if its a COMPTONIZED_PHOTON photon thats survived and change it to an injected type
                
                //replace the potantial null photon with this photon's data
                (*ph_orig)[count].p0=ph->p0;
                (*ph_orig)[count].p1=ph->p1;
                (*ph_orig)[count].p2=ph->p2;
                (*ph_orig)[count].p3=ph->p3;
                (*ph_orig)[count].comv_p0=ph->comv_p0;
                (*ph_orig)[count].comv_p1=ph->comv_p1;
                (*ph_orig)[count].comv_p2=ph->comv_p2;
                (*ph_orig)[count].comv_p3=ph->comv_p3;
                (*ph_orig)[count].r0= ph->r0;
                (*ph_orig)[count].r1=ph->r1 ;
                (*ph_orig)[count].r2=ph->r2;
                (*ph_orig)[count].s0=ph->s0;
                (*ph_orig)[count].s1=ph->s1;
                (*ph_orig)[count].s2=ph->s2;
                (*ph_orig)[count].s3=ph->s3;
                (*ph_orig)[count].num_scatt=ph->num_scatt;
                (*ph_orig)[count].weight=ph->weight;
                (*ph_orig)[count].nearest_block_index=ph->nearest_block_index;
                (*ph_orig)[count].type=ph->type;
                
                //increment count
                count+=1;
                
                if ((ph->type == COMPTONIZED_PHOTON) || (ph->type == UNABSORBED_CS_PHOTON) )
                {
                    //if the photon is a COMPTONIZED_PHOTON phton (scattered synch photon from the current frame) or a UNABSORBED_CS_PHOTON photon (scattered synch photon) from an old frame
                    //count how many of these there are
                    *scatt_cyclosynch_num_ph+=1;
                }
                
            }
             */
        }
        /*
        else
        {
            //see if the photon was a previous INJECTED_PHOTON photon absorbed that we still have to account for in the array
            if ((ph->p0 < 0) )
            {
                //replace the potantial null photon with this photon's data
                (*ph_orig)[count].p0=ph->p0;
                (*ph_orig)[count].p1=ph->p1;
                (*ph_orig)[count].p2=ph->p2;
                (*ph_orig)[count].p3=ph->p3;
                (*ph_orig)[count].comv_p0=ph->comv_p0;
                (*ph_orig)[count].comv_p1=ph->comv_p1;
                (*ph_orig)[count].comv_p2=ph->comv_p2;
                (*ph_orig)[count].comv_p3=ph->comv_p3;
                (*ph_orig)[count].r0= ph->r0;
                (*ph_orig)[count].r1=ph->r1 ;
                (*ph_orig)[count].r2=ph->r2;
                (*ph_orig)[count].s0=ph->s0;
                (*ph_orig)[count].s1=ph->s1;
                (*ph_orig)[count].s2=ph->s2;
                (*ph_orig)[count].s3=ph->s3;
                (*ph_orig)[count].num_scatt=ph->num_scatt;
                (*ph_orig)[count].weight=ph->weight;
                (*ph_orig)[count].nearest_block_index=ph->nearest_block_index;
                (*ph_orig)[count].type=ph->type;
                
                //increment count
                count+=1;
            }
        }
            */
        //fprintf(fPtr, "photon %d has energy %e and weight %e with FLASH grid number %d\n", i, ph->p0*C_LIGHT/1.6e-9, ph->weight, ph->nearest_block_index);
    }
    //fprintf(fPtr, "In phAbsCyclosynch func: abs_ph_count: %d synch_ph_count: %d scatt_cyclosynch_num_ph: %d\n", abs_ph_count, synch_ph_count, *scatt_cyclosynch_num_ph);
    *num_abs_ph=abs_ph_count; //+synch_ph_count; dont need this
    
    //fprintf(fPtr, "In phAbsCyclosynch func: count before_loop= %d\n", count);
    /*
    while (count<*num_ph)
    {
        //overwrite the last few photons to make sure that they are null photons
        (*ph_orig)[count].weight=0;
        (*ph_orig)[count].nearest_block_index=-1;
        //fprintf(fPtr, "photon %d has frequency %e and weight %e with FLASH grid number %d\n", count, (*ph_orig)[count].comv_p0*C_LIGHT/PL_CONST, (*ph_orig)[count].weight, (*ph_orig)[count].nearest_block_index);
        //fflush(fPtr);
        
        count+=1;
    }
     */ //not necessary since all photons are set to NULL with the photonList struct
    //fprintf(fPtr, "In phAbsCyclosynch func: count after loop= %d\n", count);

    return abs_count;
}


