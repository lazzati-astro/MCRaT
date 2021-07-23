void modifyRikenHydroName(char flash_file[STR_BUFFER], char prefix[STR_BUFFER], int frame);

void readHydro2D(char hydro_prefix[STR_BUFFER], int frame, double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
                     double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, FILE *fPtr);

void read_hydro(char hydro_prefix[STR_BUFFER], int frame, double r_inj, double **x, double **y, double **z, double **szx, double **szy, double **r,\
    double **theta, double **phi, double **velx, double **vely, double **velz, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number,  int ph_inj_switch, double min_r, double max_r, double fps, FILE *fPtr);

void photonInjection3D( struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, int array_length, double fps, double theta_min, double theta_max,\
double *x, double *y, double *z, double *szx, double *szy, double *r, double *theta, double *phi, double *temps, double *vx, double *vy, double *vz, gsl_rng * rand, FILE *fPtr);

int *getIndexesForRadialRemapping(char hydro_prefix[STR_BUFFER]);
