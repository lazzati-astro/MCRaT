//define on and off for switches
#define ON  1
#define OFF 0

//define switches for MCRaT to continue a simulation 'c' or initalize a simulation
#define INITALIZE 'i'
#define CONTINUE  'c'

//define the codes for the different types of hydro simulations that we can use
#define FLASH           0
#define PLUTO_CHOMBO    1
#define PLUTO           2 //have separate compiler directive for file type
#define RIKEN           3 //need to replace this when done testing, just keeping it in the background for now and minimize its effect

//define PLUTO file types here
#define FILE_DBL         0
#define FILE_FLT         1
#define FILE_DBL_H5      2
#define FILE_FLT_H5      3
#define FILE_VTK         4

//define types of simulations that can be run
#define SCIENCE                         0
#define CYLINDRICAL_OUTFLOW             1
#define SPHERICAL_OUTFLOW               2
#define STRUCTURED_SPHERICAL_OUTFLOW    3

//define the geometries that we can handle
#define CARTESIAN   0 //2D, 3D
#define SPHERICAL   1 //2D, 3D
#define CYLINDRICAL 2 //2D
#define POLAR       3 //only in 3D, technically Cylindrical coordinates but PLUTO calls this POLAR

//define the types of things that we can assume for the thermal synchrotron emission and how we calculate the B field
#define INTERNAL_E  0
#define TOTAL_E     1
#define SIMULATION  2 //need to add this to statement below to take care of defaults

//define photon types
#define INJECTED_PHOTON 'i'
#define COMPTONIZED_PHOTON 'k'
#define CS_POOL_PHOTON 'p'  
#define UNABSORBED_CS_PHOTON 'c'
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
    char type; //was the photon injected as blackbody or wien, 'i', or was it emitted as cyclo-synchrotron or was it a cyclo-synchrotron photon that was compton scattered
    double p0; //E/c, 4 momentum is in lab frame
    double p1; // p_x
    double p2; //p_y
    double p3; //p_z
    double comv_p0; //E/c, 4 momentum is in comoving frame
    double comv_p1; // p_x
    double comv_p2; //p_y
    double comv_p3; //p_z
    double r0; //x in MCRaT coordinates these get saved
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

struct hydro_dataframe
{
    /*
     Coordinate System  |   Coordinate unit vector order (r0,r1,r2)/(v0,v1,v2)
     3D Cartesian       |       x, y, z
     3D Spherical       |       r, theta, phi
     3D Polar           |       r, phi, z
     2D Cartesian       |       x, z
     2D Cylindrical     |       r, z //in PLUTO its possible to save 3D vectors, dont support this for now assume only 2D vectors
     2D Spherical       |       r, theta
     */
    int num_elements; //number of elements in each array
    double *r0; //coodinates in hydro coodinate system that user provides,
    double *r1;
    double *r2;
    double *r0_size;//size of fluid elements
    double *r1_size;
    double *r2_size;
    double *r;
    double *theta;
    double *v0; //velocity in hydro coordinate system
    double *v1;
    double *v2;
    double *dens;
    double *dens_lab;
    double *pres;
    double *temp;
    double *gamma;
    double *B0; //magentic field in hydro coordinate system
    double *B1;
    double *B2;
    
    //also hold global simulation information
    double r0_domain[2]; //holds min and max values of r0 coordinates for hydro sim
    double r1_domain[2];
    double r2_domain[2];
    double fps; //frames per second of the simulation
    int scatt_frame_number;
    int inj_frame_number;
    int last_frame;
    int increment_inj_frame; //the change in between each injection frame which may change if different number of fps is used in each portion of hydro code
    int increment_scatt_frame; //same as above expect for frames that photons are acattered in
}; // structure to hold all information for a given hydro simulation

#include "mcrat_input.h"

#include "mclib_3d.h"
#include "mclib_pluto.h"
#include "mc_cyclosynch.h"

//take care of default PLUTO file types here
#ifdef PLUTO
    #ifndef PLUTO_FILETYPE
        #define PLUTO_FILETYPE  FILE_DBL
    #endif
#endif

//take care of synchrotron defaults here
#ifdef CYCLOSYNCHROTRON_SWITCH

    //if the percentage of max photon that will be used to create the energy bins isnt defined, define it to be 10%, also applies to emiting synch photons
    #ifndef CYCLOSYNCHROTRON_REBIN_E_PERC
        #define CYCLOSYNCHROTRON_REBIN_E_PERC 0.1
    #endif

    //if the angle bins that the rebinned synch photons isnt defined use 0.5 degree increments by default
    #ifndef CYCLOSYNCHROTRON_REBIN_ANG
        #define CYCLOSYNCHROTRON_REBIN_ANG 0.5
    #endif

    //if the user hasnt defined anything for how to calculate the B field, assume that they want it calculated from the total energy
    #ifndef B_FIELD_CALC
        #define B_FIELD_CALC TOTAL_E
    #endif
    //it is defined therefore see if EPSILON_B has been set and B_FIELD_CALC != SIMULATION
    #if B_FIELD_CALC == TOTAL_E || B_FIELD_CALC == INTERNAL_E
        //see if epsilon_b has been set
        #ifndef EPSILON_B
            //if not set it to be 0.5 by default
            #define EPSILON_B 0.5
        #endif
    #endif
#else
    //if its not defined set it to be off by default
    #define CYCLOSYNCHROTRON_SWITCH OFF

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
#error Need to define hydro simulation type in mcrat_input.h file using SIM_SWITCH
#endif

#ifndef DIMENSIONS
#error Need to define hydro simulation dimensions in mcrat_input.h file using DIMENSIONS (should be set to two for this version of MCRaT)
#endif

#ifndef GEOMETRY
#error Need to define hydro simulation geometry in mcrat_input.h file using GEOMETRY
#endif

#ifndef HYDRO_L_SCALE
#error Need to define hydro simulation length scaling in mcrat_input.h file using HYDRO_L_SCALE
#endif

#ifndef HYDRO_P_SCALE
#error Need to define hydro simulation pressure scaling in mcrat_input.h file using HYDRO_P_SCALE
#endif

#ifndef HYDRO_D_SCALE
#error Need to define hydro simulation density scaling in mcrat_input.h file using HYDRO_D_SCALE
#endif

#if B_FIELD_CALC==SIMULATION
    #ifndef HYDRO_B_SCALE
    #error Need to define hydro simulation magnetic field scaling in mcrat_input.h file using HYDRO_B_SCALE
    #endif
#endif

#ifndef MCPAR
#error Need to define name of MCRaT parameter file in mcrat_input.h file using MCPAR (it is typically called mc.par, see e.g. the MCRaT manual)
#endif


void printPhotons(struct photon *ph, int num_ph, int num_ph_abs, int num_cyclosynch_ph_emit, int num_null_ph, int scatt_cyclosynch_num_ph, int frame,int frame_inj, int frame_last, char dir[STR_BUFFER], int angle_rank, FILE *fPtr );

void readMcPar(char file[STR_BUFFER], struct hydro_dataframe *hydro_data, double *theta_jmin, double *theta_j, double *n_theta_j, double **inj_radius, int **frm0, int **frm2, int *min_photons, int *max_photons, char *spect, char *restart);

void readAndDecimate(char flash_file[STR_BUFFER], struct hydro_dataframe *hydro_data, double r_inj, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);
 
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

int saveCheckpoint(char dir[STR_BUFFER], int frame,  int frame2, int scatt_frame, int ph_num,double time_now, struct photon *ph , int last_frame, int angle_rank, int angle_size);

int readCheckpoint(char dir[STR_BUFFER], struct photon **ph,  int *frame2, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time, int angle_rank, int *angle_size );

void dirFileMerge(char dir[STR_BUFFER], int start_frame, int last_frame, int numprocs,  int angle_id, FILE *fPtr);

void cylindricalPrep(double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array);

void sphericalPrep(double *r,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr);

void structuredFireballPrep(double *r, double *theta,  double *x, double *y, double *gamma, double *vx, double *vy, double *dens, double *dens_lab, double *pres, double *temp, int num_array, FILE *fPtr);

void modifyFlashName(char flash_file[STR_BUFFER], char prefix[STR_BUFFER], int frame);


void readHydro2D(char hydro_prefix[STR_BUFFER], int frame, double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
                     double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, FILE *fPtr);

int getOrigNumProcesses(int *counted_cont_procs,  int **proc_array, char dir[STR_BUFFER], int angle_rank,  int angle_procs, int last_frame);

void mcratCoordinateToHydroCoordinate(double *ph_hydro_coord, double mcrat_r0, double mcrat_r1, double mcrat_r2);

void hydroCoordinateToSpherical(double *r, double *theta, double r0, double r1, double r2);

void fillHydroCoordinateToSpherical(struct hydro_dataframe *hydro_data);

void hydroVectorToCartesian(double *cartesian_vector_3d, double v0, double v1, double v2, double x0, double x1, double x2);

void hydroDataFrameInitialize(struct hydro_dataframe *hydro_data);

void freeHydroDataFrame(struct hydro_dataframe *hydro_data);
