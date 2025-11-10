 
void photonInjection(struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, double theta_min, double theta_max, struct hydro_dataframe *hydro_data, gsl_rng * rand, FILE *fPtr);

void lorentzBoost(double *boost, double *p_ph, double *result, char object,  FILE *fPtr);

double *zeroNorm(double *p_ph);

int findContainingHydroCell( struct photon *ph, int num_ph, struct hydro_dataframe *hydro_data, int find_nearest_block_switch, FILE *fPtr);

void calcMeanFreePath(struct photon *ph, int num_ph, double *all_time_steps, int *sorted_indexes, struct hydro_dataframe *hydro_data, gsl_rng * rand, FILE *fPtr);

void reverseSortIndexes(void *sorted_indexes, int num_elements, size_t element_size, void *context_array);

int compare1(void *ar, const void *a, const void *b);

int compare2( const void *a, const void *b, void *ar);
                                   
int interpolatePropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *n_dens_lab, double *n_vx, double *n_vy, double *n_vz, double *n_temp, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr);
    
void updatePhotonPosition(struct photon *ph, int num_ph, double t, FILE *fPtr);

double photonEvent(struct photon *ph, int num_ph, double dt_max, double *all_time_steps, int *sorted_indexes, struct hydro_dataframe *hydro_data, int *scattered_ph_index, int *frame_scatt_cnt, int *frame_abs_cnt,  gsl_rng * rand, FILE *fPtr);

double averagePhotonEnergy(struct photon *ph, int num_ph);

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg, double *r_avg, FILE *fPtr  );

void phMinMax(struct photon *ph, int ph_num, double *min, double *max, double *min_theta, double *max_theta, FILE *fPtr);

