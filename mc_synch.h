/*
This header file is for the different functions for emitting and absorbing synchrotron photons
*/

double calcCyclotronFreq(double magnetic_field);

double calcDimlessTheta(double temp);

double calcB(double el_dens, double temp, double epsilon_b);

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

double synCrossSection(double el_dens, double T, double nu_ph, double p_el, double epsilon_b);

double calcSynchRLimits(int frame_scatt, int frame_inj, double fps,  double r_inj, char *min_or_max);

int photonEmitSynch(struct photon **ph_orig, int *num_ph, int *num_null_ph, double **all_time_steps, int **sorted_indexes, double r_inj, double ph_weight, int maximum_photons, int array_length, double fps, double theta_min, double theta_max , int frame_scatt, int frame_inj, double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temp, double *dens, double *vx, double *vy,  double epsilon_b, gsl_rng *rand, int riken_switch, int inject_single_switch, int scatt_ph_index, double nu_c_scatt, FILE *fPtr);

int phAbsSynch(struct photon **ph_orig, int *num_ph, int *num_abs_ph, double epsilon_b, double *temp, double *dens);
