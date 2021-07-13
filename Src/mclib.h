//define on and off for switches
#define ON  1
#define OFF 0

//define the codes for the different types of hydro simulations that we can use
#define FLASH           0
#define PLUTO_CHOMBO    1
#define PLUTO           2 //have separate compiler directive for file type
#define RIKEN           3 //need to replace this when done testing, just keeping it in the background for now and minimize its effect

//define types of simulations that can be run
#define SCIENCE                         0
#define CYLINDRICAL_OUTFLOW             1
#define SPHERICAL_OUTFLOW               2
#define STRUCTURED_SPHERICAL_OUTFLOW    3

//define the geometries that we can handle
#define CARTESIAN   0
#define SPHERICAL   1
#define CYLINDRICAL 2
#define POLAR       3

//define the types of things that we can assume for the thermal synchrotron emission and how we calculate the B field
#define INTERNAL_E  0
#define TOTAL_E     1

//define photon types
#define INJECTED_PHOTON 'i'
#define COMPTONIZED_PHOTON 'c'
#define SYNCHROTRON_POOL_PHOTON 'p'
#define OLD_COMPTONIZED_PHOTON 's'
#define REBINNED_PHOTON 'r'


extern const double C_LIGHT;
extern const double A_RAD;
extern const double PL_CONST;
extern const double K_B;
extern const double M_P;
extern const double THOM_X_SECT;
extern const double M_EL;
extern const double FINE_STRUCT;
extern const double CHARGE_EL;
extern const double R_EL;


#define STR_BUFFER 2000

struct photon
{
    char type; //was the photon injected as blackbody or wien, 'i', or was it emitted as synchrotron 's', or 'c' is it was a synchrotron photon that was compton scattered
    double p0; //E/c, 4 momentum is in lab frame
    double p1; // p_x
    double p2; //p_y
    double p3; //p_z
    double comv_p0; //E/c, 4 momentum is in comoving frame
    double comv_p1; // p_x
    double comv_p2; //p_y
    double comv_p3; //p_z
    double r0; //x in MCRaT coordinates
    double r1; //y
    double r2; //z
    double s0; //stokes I always 1
    double s1; //stokes Q/I +1 is in positive y_tilde coordinate (z_hat X photon 4 momentum)
    double s2; //stokes U/I
    double s3; //Stokes V/I
    double num_scatt;
    double weight; //each photon should have equal weight, sp this shouldnt matter, weight in mc.par file but across injections can have varying weights
    int nearest_block_index; //index that  allows for extraction of information of the hydro grid block closest to the photon
} ; //structure to hold photon information

#include "mcrat_input.h"

#include "mclib_3d.h"
#include "mclib_pluto.h"
#include "mc_synch.h"

//take care of synchrotron defaults here
#ifdef SYNCHROTRON_SWITCH

    //if the percentage of max photon that will be used to create the energy bins isnt defined, define it to be 10%, also applies to emiting synch photons
    #ifndef SYNCHROTRON_REBIN_E_PERC
        #define SYNCHROTRON_REBIN_E_PERC 0.1
    #endif

    //if the angle bins that the rebinned synch photons isnt defined use 0.5 degree increments by default
    #ifndef SYNCHROTRON_REBIN_ANG
        #define SYNCHROTRON_REBIN_ANG 0.5
    #endif

    //if the user hasnt defined anything for how to calculate the B field, assume that they want it calculated from the total energy
    #ifndef B_FIELD_CALC
        #define B_FIELD_CALC TOTAL_E
    #endif
    //it is defined therefore see if TOTAL_E is what has been set
    #if B_FIELD_CALC == TOTAL_E
        //see if epsilon_b has been set
        #ifndef EPSILON_B
            //if not set it to be 0.5 by default
            #define EPSILON_B 0.5
        #endif
    #endif
#else
    //if its not defined set it to be off by default
    #define SYNCHROTRON_SWITCH OFF

#endif

//take care of stokes switch and comv_switch and save_type defaults too
#ifndef STOKES_SWITCH
    #define STOKES_SWITCH   OFF
#endif

#ifndef COMV_SWITCH
    #define COMV_SWITCH   OFF
#endif

#ifndef SAVE_TYPE
    #define SAVE_TYPE   OFF
#endif

//throw errors during compilation if other switches are not defined
#ifndef SIM_SWITCH
#error Need to define simulation type in mcrat_input.h file using SIM_SWITCH
#endif

#ifndef DIMENSIONS
#error Need to define simulation dimensions in mcrat_input.h file using DIMENSIONS (should be set to two for this version of MCRaT)
#endif

#ifndef GEOMETRY
#error Need to define simulation geometry in mcrat_input.h file using GEOMETRY
#endif

#ifndef HYDRO_L_SCALE
#error Need to define simulation length scaling in mcrat_input.h file using HYDRO_L_SCALE
#endif

#ifndef MCPAR
#error Need to define name of MCRaT parameter file in mcrat_input.h file using MCPAR (it is typically called mc.par, see e.g. the MCRaT manual)
#endif


void printPhotons(struct photon *ph, int num_ph, int num_ph_abs, int num_ph_emit, int num_null_ph, int scatt_synch_num_ph, int frame,int frame_inj, int frame_last, char dir[200], int angle_rank, FILE *fPtr );

void readMcPar(char file[200], double *fluid_domain_x, double *fluid_domain_y, double *fps, double *theta_jmin, double *theta_j, double *d_theta_j, double *inj_radius_small, double *inj_radius_large, int *frm0_small, int *frm0_large,\
int *last_frm, int *frm2_small,int *frm2_large, double *ph_weight_small,double *ph_weight_large,int *min_photons, int *max_photons, char *spect, char *restart);

void readAndDecimate(char flash_file[200], double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
 double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);
 
 void photonInjection( struct photon **ph, int *ph_num, double r_inj, double ph_weight, int min_photons, int max_photons, char spect, int array_length, double fps, double theta_min, double theta_max,\
double *x, double *y, double *szx, double *szy, double *r, double *theta, double *temps, double *vx, double *vy, gsl_rng * rand, FILE *fPtr);

void lorentzBoost(double *boost, double *p_ph, double *result, char object,  FILE *fPtr);

double *zeroNorm(double *p_ph);

int findNearestBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z);

int findContainingBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z, double *szx, double *szy, int old_block_index, int find_block_switch, FILE *fPtr);

int checkInBlock(int block_index, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z, double *szx, double *szy);

int findNearestPropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double hydro_domain_x, double hydro_domain_y, double epsilon_b, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *all_time_steps,  int *sorted_indexes, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr);
                                   
int compare (void *ar, const void *a, const void *b);

int compare2 ( const void *a, const void *b, void *ar);
                                   
int interpolatePropertiesAndMinMFP( struct photon *ph, int num_ph, int array_num, double *time_step, double *x, double  *y, double *z, double *szx, double *szy, double *velx,  double *vely, double *velz, double *dens_lab,\
                                   double *temp, double *n_dens_lab, double *n_vx, double *n_vy, double *n_vz, double *n_temp, gsl_rng * rand, int find_nearest_block_switch, FILE *fPtr);
    
void updatePhotonPosition(struct photon *ph, int num_ph, double t, FILE *fPtr);

void mullerMatrixRotation(double theta, double *s, FILE *fPtr);

void findXY(double *v_ph, double *vector, double *x, double *y);

void stokesRotation(double *v, double *v_ph, double *v_ph_boosted, double *s, FILE *fPtr);

double photonEvent(struct photon *ph, int num_ph, double dt_max, double *all_time_steps, int *sorted_indexes, double *all_flash_vx, double *all_flash_vy, double *all_flash_vz, double *all_fluid_temp, int *scattered_ph_index, int *frame_scatt_cnt, int *frame_abs_cnt, gsl_rng * rand, FILE *fPtr);

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr);

int singleScatter(double *el_comov, double *ph_comov, double *s, gsl_rng * rand, FILE *fPtr);

int comptonScatter(double *theta, double *phi, gsl_rng * rand, FILE *fPtr);

int kleinNishinaScatter(double *theta, double *phi, double p0, double q, double u, gsl_rng * rand, FILE *fPtr);

double averagePhotonEnergy(struct photon *ph, int num_ph);

void phScattStats(struct photon *ph, int ph_num, int *max, int *min, double *avg, double *r_avg, FILE *fPtr  );

int saveCheckpoint(char dir[200], int frame,  int frame2, int scatt_frame, int ph_num,double time_now, struct photon *ph , int last_frame, int angle_rank, int angle_size);

int readCheckpoint(char dir[200], struct photon **ph,  int *frame2, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time, int angle_rank, int *angle_size );

void dirFileMerge(char dir[200], int start_frame, int last_frame, int numprocs,  int angle_id, FILE *fPtr);

void cylindricalPrep(double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array);

void sphericalPrep(double *r,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr);

void structuredFireballPrep(double *r, double *theta,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr);

void modifyFlashName(char flash_file[200], char prefix[200], int frame);


void readHydro2D(char hydro_prefix[200], int frame, double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
                     double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, FILE *fPtr);

int getOrigNumProcesses(int *counted_cont_procs,  int **proc_array, char dir[200], int angle_rank,  int angle_procs, int last_frame);
