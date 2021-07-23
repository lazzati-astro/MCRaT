//
//  mcrat_io.h
//  contains functions that handle IO for MCRaT
//  includes wrapper functions for reading in hydro files and functios for MCRaT to save its own data
//  Created by Tyler Parsotan on 7/23/21.
//

int getOrigNumProcesses(int *counted_cont_procs,  int **proc_array, char dir[STR_BUFFER], int angle_rank,  int angle_procs, int last_frame);

void printPhotons(struct photon *ph, int num_ph, int num_ph_abs, int num_cyclosynch_ph_emit, int num_null_ph, int scatt_cyclosynch_num_ph, int frame,int frame_inj, int frame_last, char dir[STR_BUFFER], int angle_rank, FILE *fPtr );

int saveCheckpoint(char dir[STR_BUFFER], int frame,  int frame2, int scatt_frame, int ph_num,double time_now, struct photon *ph , int last_frame, int angle_rank, int angle_size);

int readCheckpoint(char dir[STR_BUFFER], struct photon **ph,  int *frame2, int *framestart, int *scatt_framestart, int *ph_num, char *restart, double *time, int angle_rank, int *angle_size );

void readMcPar(struct hydro_dataframe *hydro_data, double *theta_jmin, double *theta_j, double *n_theta_j, double **inj_radius, int **frm0, int **frm2, int *min_photons, int *max_photons, char *spect, char *restart);

void dirFileMerge(char dir[STR_BUFFER], int start_frame, int last_frame, int numprocs,  int angle_id, FILE *fPtr);

void hydroDataFrameInitialize(struct hydro_dataframe *hydro_data);

void freeHydroDataFrame(struct hydro_dataframe *hydro_data);

int getHydroData(struct hydro_dataframe *hydro_data, int frame, double inj_radius, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);
