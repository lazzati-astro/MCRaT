extern const double C_LIGHT;
extern const double A_RAD;
extern const double PL_CONST;
extern const double K_B;

#define STR_BUFFER 2000
struct photon
{
    double p0; //E/c, 4 momentum is in lab frame
    double p1; // p_x
    double p2; //p_y
    double p3; //p_z
    double r0; //x in MCRaT coordinates
    double r1; //y
    double r2; //z
    double s0; //stokes I always 1
    double s1; //stokes Q/I 
    double s2; //stokes U/I
    double s3; //Stokes V/I
    double num_scatt;
    double weight; //each photon should have equal weight, sp this shouldnt matter, weight in mc.par file but across injections can have varying weights
    int nearest_block_index; //index that  allows for extraction of information of the hydro grid block closest to the photon
} ; //structure to hold photon information


void printPhotons(struct photon *ph, int num_ph, int frame, int  frame_inj,char dir[200], int angle_rank, int stokes_switch, FILE *fPtr );

void readMcPar(char file[200], double *fluid_domain_x, double *fluid_domain_y, double *fps, double *theta_jmin, double *theta_j, double *d_theta_j, double *inj_radius_small, double *inj_radius_large, int *frm0_small, int *frm0_large,\
int *last_frm, int *frm2_small,int *frm2_large, double *ph_weight_small,double *ph_weight_large,int *min_photons, int *max_photons, char *spect, char *restart,  int *dim_switch);

void readAndDecimate(char flash_file[200], double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
 double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);
 
 void photonInjection( struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, int array_length, double fps, double theta_min, double theta_max,\
double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temps, double *vx, double *vy, gsl_rng * rand,  int riken_switch, FILE *fPtr);

void lorentzBoost(double *boost, double *p_ph, double *result, char object,  FILE *fPtr);

double *zeroNorm(double *p_ph);

int findNearestBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z,  int dim_switch_3d);

int findContainingBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z, double *szx, double *szy, int old_block_index, int find_block_switch, int dim_switch_3d, int riken_switch, FILE *fPtr);

int checkInBlock(int block_index, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z, double *szx, double *szy, int dim_switch_3d, int riken_switch);

int findNearestPropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double hydro_domain_x, double hydro_domain_y,  double *time_step, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *all_time_steps, int *sorted_indexes, gsl_rng * rand, int dim_switch_3d, int find_nearest_block_switch, int riken_switch, FILE *fPtr);
                                   
int compare (void *ar, const void *a, const void *b);

int compare2 ( const void *a, const void *b, void *ar);
                                   
int interpolatePropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *n_dens_lab, double *n_vx, double *n_vy, double *n_vz, double *n_temp, gsl_rng * rand, int dim_switch_3d, int find_nearest_block_switch, int riken_switch, FILE *fPtr);
    
void updatePhotonPosition(struct photon *ph, int num_ph, double t, FILE *fPtr);

void mullerMatrixRotation(double theta, double *s, FILE *fPtr);

void findXY(double *v_ph, double *vector, double *x, double *y);

void stokesRotation(double *v, double *v_ph, double *v_ph_boosted, double *s, FILE *fPtr);

double photonScatter(struct photon *ph, int num_ph, double dt_max, double *all_time_steps, int *sorted_indexes, double *all_flash_vx, double *all_flash_vy, double *all_flash_vz, double *all_fluid_temp, int *scattered_ph_index, int *frame_scatt_cnt, gsl_rng * rand, int dim_switch_3d, int stokes_switch, FILE *fPtr);

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr);

int singleScatter(double *el_comov, double *ph_comov, double *s, gsl_rng * rand, int stokes_switch, FILE *fPtr);

int comptonScatter(double *theta, double *phi, gsl_rng * rand, FILE *fPtr);

int kleinNishinaScatter(double *theta, double *phi, double p0, double q, double u, gsl_rng * rand, int stokes_switch, FILE *fPtr);

double averagePhotonEnergy(struct photon *ph, int num_ph);

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg, double *r_avg  );

int saveCheckpoint(char dir[200], int frame,  int frame2, int scatt_frame, int ph_num,double time_now, struct photon *ph , int last_frame, int angle_rank, int angle_size);

void readCheckpoint(char dir[200], struct photon **ph,  int *frame2, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time, int angle_rank, int *angle_size, int dim_switch, int riken_switch );

void dirFileMerge(char dir[200], int start_frame, int last_frame, int numprocs,  int angle_id, int dim_switch, int riken_switch, FILE *fPtr);

void cylindricalPrep(double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array);

void sphericalPrep(double *r,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr);

void structuredFireballPrep(double *r, double *theta,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr);

void modifyFlashName(char flash_file[200], char prefix[200], int frame, int dim_switch);


void readHydro2D(char hydro_prefix[200], int frame, double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
                     double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, FILE *fPtr);

int getOrigNumProcesses(int *counted_cont_procs,  int **proc_array, char dir[200], int angle_rank,  int angle_procs, int last_frame, int dim_switch, int riken_switch);
