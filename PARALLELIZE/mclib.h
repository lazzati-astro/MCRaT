struct photon
{
    double p0; //E/c, 4 momentum is in lab frame
    double p1; // p_x
    double p2; //p_y
    double p3; //p_z
    double r0; //x in MCRaT coordinates
    double r1; //y
    double r2; //z
    double num_scatt;
    //double weight; //each photon should have equal weight, sp this shouldnt matter, weight in mc.par file
} ; //structure to hold photon information


void printPhotons(struct photon *ph, int num_ph, int frame,char dir[200] );

void readMcPar(char file[200], double *fps, double *theta_jmin, double *theta_j, double *inj_radius, int *frm0,int *last_frm, int *frm2, int *photon_num, double *ph_weight, char *spect, char *restart);

void readAndDecimate(char flash_file[200], double r_inj, double **x, double **y, double **szx, double **szy, double **r,\
 double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number);
 
 void photonInjection( struct photon **ph, int *ph_num, double r_inj, double ph_weight, char spect, int array_length, double fps, double theta_min, double theta_max,\
double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temps, double *vx, double *vy, gsl_rng * rand);

void lorentzBoost(double *boost, double *p_ph, double *result, char object);

double *zeroNorm(double *p_ph);

int findNearestPropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *velx,  double *vely, double *dens_lab,\
    double *temp, double *n_dens_lab, double *n_vx, double *n_vy,double *n_temp, gsl_rng * rand);
    
void updatePhotonPosition(struct photon *ph, int num_ph, double t);

void photonScatter(struct photon *ph, double flash_vx, double flash_vy, double fluid_temp, gsl_rng * rand);

void singleElectron(double *el_p, double temp, struct photon *ph, gsl_rng * rand);

void singleComptonScatter(double *el_comov, double *ph_comov, gsl_rng * rand);

double averagePhotonEnergy(struct photon *ph, int num_ph);

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg );

void saveCheckpoint(char dir[200], int frame, int scatt_frame, int ph_num,double time_now, struct photon *ph , int last_frame);

void readCheckpoint(char dir[200], struct photon **ph, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time );
