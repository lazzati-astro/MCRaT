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
#define CUSTOM_OUTFLOW                  4

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
    int nearest_block_index; //index that  allows for extraction of information of the hydro grid block that the photon si located within
} ; //structure to hold photon information

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
}; // structure to hold all information for a given hydro simulation

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
#include <sys/types.h>
#include <sys/stat.h>
#include <glob.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <time.h>
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

//include all other header files in mcrat
#include "mcrat_input.h"

#include "mclib.h"
#include "mclib_riken.h"
#include "mclib_pluto.h"
#include "mclib_flash.h"
#include "mc_cyclosynch.h"
#include "geometry.h"
#include "mcrat_scattering.h"
#include "mcrat_io.h"
#include "custom_outflow.h"
#include "custom_spectrum.h"

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

