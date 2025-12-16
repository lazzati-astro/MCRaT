//
//  geometry.h
//  
//
//  Created by Tyler Parsotan on 7/23/21.
//


void mcratCoordinateToHydroCoordinate(double *ph_hydro_coord, double mcrat_r0, double mcrat_r1, double mcrat_r2);

void hydroCoordinateToSpherical(double *r, double *theta, double r0, double r1, double r2);

void hydroCoordinateToMcratCoordinate(double *hydro_mcrat_coord, double hydro_r0, double hydro_r1, double hydro_r2);

void fillHydroCoordinateToSpherical(struct hydro_dataframe *hydro_data);

double vectorMagnitude(double v0, double v1, double v2);

void hydroVectorToCartesian(double *cartesian_vector_3d, double v0, double v1, double v2, double x0, double x1, double x2);

double hydroElementVolume(struct hydro_dataframe *hydro_data, int index);

int findNearestBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z);

int findContainingBlock(double ph_hydro_r0, double ph_hydro_r1, double ph_hydro_r2, struct hydro_dataframe *hydro_data, FILE *fPtr);

int checkInBlock(double ph_hydro_r0, double ph_hydro_r1, double ph_hydro_r2, struct hydro_dataframe *hydro_data, int block_index);

int findContainingBlock_grid(double ph_hydro_r0, double ph_hydro_r1, double ph_hydro_r2, struct hydro_dataframe *hydro_data, FILE *fPtr);

void freeSpatialGrid(struct SpatialGrid *g);

struct SpatialGrid *buildSpatialGrid(struct hydro_dataframe *hydro_data, FILE *fPtr);
