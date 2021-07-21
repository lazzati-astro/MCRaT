/*
This header file is for the different functions for emitting and absorbing synchrotron photons
*/

double calcCyclotronFreq(double magnetic_field);

double calcEB(double magnetic_field);

double calcBoundaryE(double magnetic_field, double temp);

double calcDimlessTheta(double temp);

double calcB(double el_dens, double temp);

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

int rebinCyclosynchCompPhotons(struct photon **ph_orig, int *num_ph, int *num_null_ph, int *num_cyclosynch_ph_emit, int *scatt_cyclosynch_num_ph, double **all_time_steps, int **sorted_indexes, int max_photons, double thread_theta_min, double thread_theta_max , gsl_rng * rand, FILE *fPtr);

int rebin2dCyclosynchCompPhotons(struct photon **ph_orig, int *num_ph,  int *num_null_ph, int *num_cyclosynch_ph_emit, int *scatt_cyclosynch_num_ph, double **all_time_steps, int **sorted_indexes, int max_photons, double thread_theta_min, double thread_theta_max , gsl_rng * rand, FILE *fPtr);

int photonEmitCyclosynch(struct photon **ph_orig, int *num_ph, int *num_null_ph, double **all_time_steps, int **sorted_indexes, double r_inj, double ph_weight, int maximum_photons, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng *rand, int inject_single_switch, int scatt_ph_index, FILE *fPtr);

double phAbsCyclosynch(struct photon **ph_orig, int *num_ph, int *num_abs_ph, int *scatt_cyclosynch_num_ph, double *temp, double *dens, FILE *fPtr);

