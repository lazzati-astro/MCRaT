//
//  mcrat_input.h
//  This file serves as an input of the various types of physics that MCRaT should consider.
//  as well as the microphysical considerations that go along with those physics. Various
//  options are outlined in the documentation.
//
//  Created by Tyler Parsotan on 12/9/19.
//

#define SIMULATION_TYPE SPHERICAL_OUTFLOW
//#define FILEPATH "/Users/Tylerparsotan/Documents/16OI_TEST/"
//#define FILEROOT "rhd_jet_big_16OI_hdf5_plt_cnt_"
//#define MC_PATH "SKN_16OI_SPHERICAL_flash/"
//#define     SIM_SWITCH                  FLASH
//#define     GEOMETRY                    CARTESIAN
//#define     DIMENSIONS                  TWO

#define FILEPATH "/Users/Tylerparsotan/Documents/CHOMBO_TEST/CHOMBO_2D/"
#define FILEROOT "data."
#define MC_PATH "MCRAT_TEST/"
#define     SIM_SWITCH                  PLUTO_CHOMBO
#define     GEOMETRY                    SPHERICAL
#define     DIMENSIONS                  TWO

#define     STOKES_SWITCH               ON
#define     COMV_SWITCH                 ON
#define     HYDRO_L_SCALE               1e9
#define     HYDRO_P_SCALE               C_LIGHT*C_LIGHT
#define     HYDRO_D_SCALE               1
#define     CYCLOSYNCHROTRON_SWITCH     ON
#define     SAVE_TYPE                   ON

#define     MCPAR                   "mc.par"

/*
 Modify parameters above this comment only
 */
