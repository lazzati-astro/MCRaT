//
//  mcrat_input.h
//  This file serves as an input of the various types of physics that MCRaT should consider.
//  as well as the microphysical considerations that go along with those physics. Various
//  options are outlined in the documentation.
//
//  Created by Tyler Parsotan on 12/9/19.
//

//#define SIMULATION_TYPE CYLINDRICAL_OUTFLOW
//#define FILEPATH "/Users/Tylerparsotan/Documents/HYDRO_SIMS/16OI_TEST/"
//#define FILEROOT "rhd_jet_big_16OI_hdf5_plt_cnt_"
//#define MC_PATH "SKN_16OI_CYLINDRICAL/"
//#define     SIM_SWITCH                  FLASH
//#define     GEOMETRY                    CARTESIAN
//#define     DIMENSIONS                  TWO
//#define     HYDRO_L_SCALE               1e9


//#define SIMULATION_TYPE CYLINDRICAL_OUTFLOW
//#define FILEPATH "/Users/Tylerparsotan/Documents/HYDRO_SIMS/CHOMBO_TEST/CHOMBO_2D/"
//#define FILEROOT "data."
//#define MC_PATH "MCRAT_TEST/"
//#define     SIM_SWITCH                  PLUTO_CHOMBO
//#define     GEOMETRY                    SPHERICAL
//#define     DIMENSIONS                  TWO
//#define     HYDRO_L_SCALE               1e11

#define SIMULATION_TYPE CYLINDRICAL_OUTFLOW
#define FILEPATH "/Users/Tylerparsotan/Documents/HYDRO_SIMS/CHOMBO_TEST/CHOMBO_3D_HIGH_RES/"
#define FILEROOT "data."
#define MC_PATH "MCRAT_TEST/"
#define     SIM_SWITCH                  PLUTO_CHOMBO
#define     GEOMETRY                    SPHERICAL
#define     DIMENSIONS                  THREE
#define     HYDRO_L_SCALE               1e12

//#define SIMULATION_TYPE CYLINDRICAL_OUTFLOW
//#define FILEPATH "/Users/Tylerparsotan/Documents/HYDRO_SIMS/3D_MHD_SIMULATION_FRAMES/"
//#define FILEROOT "data."
//#define MC_PATH "MCRAT_TEST/"
//#define     SIM_SWITCH                  PLUTO
//#define     GEOMETRY                    CARTESIAN
//#define     DIMENSIONS                  THREE
////#define     B_FIELD_CALC                SIMULATION
//#define     HYDRO_L_SCALE               1e11


//#define SIMULATION_TYPE CYLINDRICAL_OUTFLOW
//#define FILEPATH "/Users/Tylerparsotan/Documents/HYDRO_SIMS/LEO_2.5D_MHD_PLUTO/BPT5/"
//#define FILEROOT "data."
//#define MC_PATH "MCRAT_TEST/"
//#define     SIM_SWITCH                  PLUTO
//#define     GEOMETRY                    CYLINDRICAL
//#define     DIMENSIONS                  TWO_POINT_FIVE
////#define     B_FIELD_CALC                SIMULATION
//#define CYCLOSYNCHROTRON_REBIN_E_PERC 0.5
//#define EPSILON_B 1.0
//#define     HYDRO_L_SCALE               1e12




#define     STOKES_SWITCH               ON
#define     COMV_SWITCH                 ON
#define     HYDRO_D_SCALE               1

#define     CYCLOSYNCHROTRON_SWITCH     OFF
#define     SAVE_TYPE                   ON

#define     MCPAR                   "mc.par"

/*
 Modify parameters above this comment only
 */
