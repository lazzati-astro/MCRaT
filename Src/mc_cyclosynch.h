/*
This header file is for the different functions for emitting and absorbing synchrotron photons
*/

/* Data structure to hold photon range information */
struct PhotonRangeInfo
{
    double p0_min;
    double p0_max;
    double log_p0_min;
    double log_p0_max;
    double theta_min;
    double theta_max;
    double phi_min;
    double phi_max;
    int synch_photon_count;
    int valid_photon_count;
};

/* Data structure for binning parameters */
struct BinningParams
{
    int num_bins;
    int num_bins_theta;
    int num_bins_phi;
    int num_avg;
    int total_bins;
};

/* Structure to hold bin statistics */
struct BinStats
{
    double weighted_r;
    double weighted_theta;
    double weighted_phi_offset;
    double weighted_stokes[4];
    double weighted_scatt_count;
    double total_weight;
    double weighted_phi_dir;
    double weighted_theta_dir;
    double weighted_energy;
#if DIMENSIONS == THREE
    double weighted_phi_pos;
#endif
};


double calcCyclotronFreq(double magnetic_field);

double calcEB(double magnetic_field);

double calcBoundaryE(double magnetic_field, double temp);

double calcDimlessTheta(double temp);

double calcB(double el_dens, double temp);

double getMagneticFieldMagnitude(struct hydro_dataframe *hydro_data, int hydro_grid_index);

double n_el_MJ(double el_dens, double dimlesstheta, double gamma);

double n_el_MB(double el_dens, double dimlesstheta, double gamma);

double Z(double nu, double nu_c, double gamma );

double Z_sec_der(double nu, double nu_c, double gamma);

double chi(double dimlesstheta, double gamma);

double gamma0(double nu, double nu_c, double dimlesstheta);

double jnu(double nu, double nu_c, double dimlesstheta, double el_dens);

double jnu_ph_spect(double nu, void *p);

double C(double nu_ph, double nu_c, double gamma_el, double p_el);

double G(double gamma_el, double p_el);

double G_prime(double gamma_el, double p_el);

double synCrossSection(double el_dens, double T, double nu_ph, double p_el);

double calcCyclosynchRLimits(int frame_scatt, int frame_inj, double fps,  double r_inj, char *min_or_max);

int rebinCyclosynchCompPhotons(struct photonList *photon_list, int *num_cyclosynch_ph_emit, int *scatt_cyclosynch_num_ph, int max_photons, double thread_theta_min, double thread_theta_max , gsl_rng * rand, FILE *fPtr);

int photonEmitCyclosynch(struct photonList *photon_list, double r_inj, double ph_weight, int maximum_photons, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng *rand, int inject_single_switch, int scatt_ph_index, FILE *fPtr);

double phAbsCyclosynch(struct photonList *photon_list, int *num_abs_ph, int *scatt_cyclosynch_num_ph, struct hydro_dataframe *hydro_data, FILE *fPtr);

static void calculate_photon_position(const struct photon *ph, double *r, double *theta, double *phi);

static int collect_photon_statistics(const struct photonList *photon_list, struct PhotonRangeInfo *info, FILE *fPtr);

static int allocate_histograms(gsl_histogram2d **h_energy_theta, gsl_histogram2d **h_energy_phi, gsl_histogram2d **h_theta_phi, int num_bins, int num_bins_theta, int num_bins_phi, double log_p0_min, double log_p0_max, double theta_min, double theta_max, double phi_min, double phi_max);

static void free_histograms(gsl_histogram2d *h_energy_theta, gsl_histogram2d *h_energy_phi, gsl_histogram2d *h_theta_phi);

static struct BinStats* allocate_bin_stats(int total_bins, int num_avg);

static void free_bin_stats(struct BinStats *stats, int total_bins);

static int calculate_bin_index(int count_x, int count_y, int count_z, int num_bins, int num_bins_theta, int num_bins_phi);
