//
//  mcrat.h
//  Global header file for all of MCRaT
//
//  Created by Tyler Parsotan on 7/23/21.
//

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

//define the hydro dimensions that we can handle
#define TWO                 0
#define TWO_POINT_FIVE      1
#define THREE               2

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
#define NULL_PHOTON 'N'

//define nonthermal functions for electron distribution
#define POWERLAW 1
#define BROKENPOWERLAW 2

//define ways that the optical depth can be calculated (eithe directly from the fluid or from the pretabulated values)
#define DIRECT 1
#define TABLE 2

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

//include all libraries needed
#include "hdf5.h"
#include "mpi.h"

#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <glob.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

//include all other header files in mcrat
#include "mcrat_input.h"

//set the nonthermal electrons to be off
//this is placed here so we can define it in prep for the photon struct
#ifndef NONTHERMAL_E_DIST
    #define NONTHERMAL_E_DIST OFF
#endif

//then include this file also in prep for the photon struct optical depth array being defined
#include "hot_x_section.h"


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
    int recalc_properties; //if this is 1 then the photon has scattered and we need to recalc the optical depths or the photon is in a new hydro cell and we need to recalc the optical depths
    double weight; //each photon should have equal weight, sp this shouldnt matter, weight in mc.par file but across injections can have varying weights
    int nearest_block_index; //index that  allows for extraction of information of the hydro grid block that the photon si located within
    double time_to_scatter; //the sampled mean free path of the photon divided by C_LIGHT
    #if NONTHERMAL_E_DIST != OFF
        double optical_depths[1+N_GAMMA]; //the optical depths that are calculated for thermal + non-thermal electrons with the nonthermal electron subgroups
        double scattering_bias[1+N_GAMMA];
    #endif
    //save the total calculated optical depth, if we only have thermal electrons this optical depth is that calcualted value otherwise it includes thermal-non-thermal electrons and the biases
    double total_optical_depth;
} ; //structure to hold photon information

struct photonList
{
    struct photon *photons;
    int *sorted_indexes
    int num_photons;
    int num_null_photons;
    int list_capacity;
};

struct hydro_dataframe
{
    /*
     Coordinate System  |   Coordinate unit vector order (r0,r1,r2)/(v0,v1,v2)
     3D Cartesian       |       x, y, z
     3D Spherical       |       r, theta, phi
     3D Polar           |       r, phi, z
     2D Cartesian       |       x, z
     2D Cylindrical     |       r, z (phi) //in PLUTO its possible to save 3D vectors  in 2.5 dims where the final unit vector is phi hat
     2D Spherical       |       r, theta, (phi)
     */
    int num_elements; //number of elements in each array
    double *r0; //coodinates in hydro coodinate system that user provides,
    double *r1;
    double *r2;
    double *r0_size;//size of fluid elements
    double *r1_size;
    double *r2_size;
    double *r;//spherical coordinates of fluid elements for photon injection, mdifying fluid value purposes
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
    #if NONTHERMAL_E_DIST != OFF
        double electron_dens_subgroup[N_GAMMA]; //the density of nonthermal electrons in each subgroup (needs to be multiplied by the number density of nonthermal electrons to get the actual number density)
        double average_dimless_theta; //the volume average dimensionless temperature for us to use in calculating scattering bias
        double *nonthermal_dens; //this is the number density of non-thermal electrons in each cell (usually defined based on B field)
    #endif

}; // structure to hold all information for a given hydro simulation


#include "mclib.h"
#include "mclib_riken.h"
#include "mclib_pluto.h"
#include "mclib_flash.h"
#include "mc_cyclosynch.h"
#include "geometry.h"
#include "mcrat_scattering.h"
#include "mcrat_io.h"
#include "analytic_outflows.h"
#include "optical_depth.h"
#include "electron.h"


//if the user doesnt specify NONTHERMAL_E_DIST set it to be off
#ifndef NONTHERMAL_E_DIST
    #define NONTHERMAL_E_DIST OFF
#endif

//if the user specifies the NONTHERMAL_E_DIST then they need to also specify
// TAU_CALCULATION = TABLE to use the tabulated cross sectional values
#if NONTHERMAL_E_DIST != OFF
    //set the default optical depth calculation to be that of the fluid properties
    #ifndef TAU_CALCULATION
        #define TAU_CALCULATION TABLE
    #endif

    //The user can specify NONTHERMAL_E_DIST and set TAU_CALCULATION = DIRECT, this is not permitted so throw an error
    #if TAU_CALCULATION == DIRECT
        #error NONTHERMAL_E_DIST cannot be set while TAU_CALCULATION = DIRECT.
    #endif
#endif

//set the default optical depth calculation to be that of the fluid properties
#ifndef TAU_CALCULATION
    #define TAU_CALCULATION DIRECT
#endif

//set default hydro v scale to be the speed of light
#define HYDRO_V_SCALE C_LIGHT

//take care of hydro scales in pressure, magnetic field, etc
#define HYDRO_P_SCALE   HYDRO_D_SCALE*HYDRO_V_SCALE*HYDRO_V_SCALE
#if B_FIELD_CALC==SIMULATION
    #define HYDRO_B_SCALE sqrt(4*M_PI*HYDRO_P_SCALE)
#endif

//take care of default PLUTO file types here
#ifdef PLUTO
    #ifndef PLUTO_FILETYPE
        #define PLUTO_FILETYPE  FILE_DBL
    #endif
#endif

#ifndef CYCLOSYNCHROTRON_SWITCH
    #define CYCLOSYNCHROTRON_SWITCH OFF
#endif

//take care of synchrotron defaults here
#ifdef CYCLOSYNCHROTRON_SWITCH

    //if the percentage of max photon that will be used to create the energy bins isnt defined, define it to be 10%, also applies to emiting synch photons
    #ifndef CYCLOSYNCHROTRON_REBIN_E_PERC
        #define CYCLOSYNCHROTRON_REBIN_E_PERC 0.1
    #endif

    //if the polar angle bins that the rebinned synch photons isnt defined use 0.5 degree increments by default
    #ifndef CYCLOSYNCHROTRON_REBIN_ANG
        #define CYCLOSYNCHROTRON_REBIN_ANG 0.5
    #endif

    //if the azimuthal angle bins that the rebinned synch photons isnt defined use 0.5 degree increments by default
    #ifndef CYCLOSYNCHROTRON_REBIN_ANG_PHI
        #define CYCLOSYNCHROTRON_REBIN_ANG_PHI 10 //0.5
    #endif

    #if NONTHERMAL_E_DIST == OFF
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
    #endif
#endif

#if NONTHERMAL_E_DIST != OFF
    //if the user is specifying one of the 2 on-thermal distributions to use need to do error checking
    #if NONTHERMAL_E_DIST == POWERLAW
        #ifndef POWERLAW_INDEX
            #error Need to define POWERLAW_INDEX in mcrat_input.h file using POWERLAW_INDEX
        #endif

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
        #ifndef POWERLAW_INDEX_1
            #error Need to define POWERLAW_INDEX_1 in mcrat_input.h file using POWERLAW_INDEX_1
        #endif

        #ifndef POWERLAW_INDEX_2
            #error Need to define POWERLAW_INDEX_2 in mcrat_input.h file using POWERLAW_INDEX_2
        #endif

        #ifndef GAMMA_BREAK
            #error Need to define GAMMA_BREAK in mcrat_input.h file using GAMMA_BREAK
        #endif

    #else
        #error Unnknown nonthermal electron distribution.
    #endif

    #if NONTHERMAL_E_DIST == POWERLAW || NONTHERMAL_E_DIST == BROKENPOWERLAW
        #ifndef GAMMA_MIN
            #error Need to define GAMMA_MIN in mcrat_input.h file using GAMMA_MIN
        #endif

        #ifndef GAMMA_MAX
            #error Need to define GAMMA_MAX in mcrat_input.h file using GAMMA_MAX
        #endif
    #endif

    //if the user hasnt defined anything for how to calculate the B field, assume that they want it calculated from the total energy
    #ifndef B_FIELD_CALC
        #error B_FIELD_CALC needs to be defined with NONTHERMAL_E_DIST. Specify B_FIELD_CALC in mcrat_input.h file using B_FIELD_CALC.
    #endif

    //it is defined therefore see if EPSILON_B has been set and B_FIELD_CALC != SIMULATION
    #if B_FIELD_CALC == TOTAL_E || B_FIELD_CALC == INTERNAL_E
    //see if epsilon_b has been set
        #ifndef EPSILON_B
            //dont assume anything here
            #error EPSILON_B needs to be defined with this B_FIELD_CALC setting. Specify EPSILON_B in mcrat_input.h file using EPSILON_B
        #endif
    #endif

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

#ifndef HYDRO_D_SCALE
#error Need to define hydro simulation density scaling in mcrat_input.h file using HYDRO_D_SCALE
#endif

#ifndef MCPAR
#error Need to define name of MCRaT parameter file in mcrat_input.h file using MCPAR (it is typically called mc.par, see e.g. the MCRaT manual)
#endif

