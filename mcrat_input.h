//
//  mcrat_input.h
//  
//
//  Created by Tyler Parsotan on 12/9/19.
//



//MODIFY C COMPILER DIRECTIVES BELOW THIS LINE
/*
#define THISRUN "Science"
#define FILEPATH "/Volumes/LACIE_RAID/Collapsars/2D/HUGE_BOXES/CONSTANT/16TI/"
#define FILEROOT "rhd_jet_big_13_hdf5_plt_cnt_"
#define MC_PATH "KN_CMC_16TI/"
*/
 #define    THISRUN                 "Science"
 #define    FILEPATH                "/Users/parsotat/Downloads/"
 #define    FILEROOT                "data."
 #define    MC_PATH                 "PLUTO_MCRAT/"
 
/*
 #define THISRUN "Science"
 #define FILEPATH "/home/physics/parsotat/16OM/"
 #define FILEROOT "rhd_jet_big_16OM_hdf5_plt_cnt_"
 #define MC_PATH "DIR_TEST/"
 
 #define THISRUN "Spherical"
 #define FILEPATH "/Volumes/DATA6TB/Collapsars/2D/HUGE_BOXES/CONSTANT/16OI/"
 //#define FILEPATH "/Users/Tylerparsotan//Documents/16OI_TEST/"
 #define FILEROOT "rhd_jet_big_16OI_hdf5_plt_cnt_"
 #define MC_PATH "TEST/"
 
 #define THISRUN "Spherical"
 //#define THISRUN "Structured Spherical"
 //#define FILEPATH "/home/physics/parsotat/16TI/"
 #define FILEPATH "/Users/Tylerparsotan/Documents/16TI/"
 #define FILEROOT "rhd_jet_big_13_hdf5_plt_cnt_"
 #define MC_PATH "KN_CMC_16TI_SPHERICAL/"
*/


#define     SIM_SWITCH              PLUTO_CHOMBO
#define     STOKES_SWITCH           ON
#define     COMV_SWITCH             OFF
#define     DIMENSIONS              2
#define     GEOMETRY                SPHERICAL
#define     HYDRO_L_SCALE           1e9

#define     MCPAR                   "mc.par"

/*
 Modify parameters above this comment only
 */
