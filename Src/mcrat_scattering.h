//
//  mcrat_scattering.h
//  
//
//  Created by Tyler Parsotan on 7/23/21.
//

void mullerMatrixRotation(double theta, double *s, FILE *fPtr);

void findXY(double *v_ph, double *vector, double *x, double *y);

double findPhi(double *x_old, double *y_old, double *x_new, double *y_new);

void stokesRotation(double *v, double *v_ph, double *v_ph_boosted, double *s, FILE *fPtr);

void stokesScatter(double *s, double *orig_s,  gsl_vector *ph_p_orig, gsl_vector *result0, double ph_p_prime, double scattered_ph_e, FILE *fPtr);

int singleScatter(double *el_comov, double *ph_comov, double *s, gsl_rng * rand, FILE *fPtr);

int comptonScatter(double *theta, double *phi, gsl_rng * rand, FILE *fPtr);

int kleinNishinaScatter(double *theta, double *phi, double p0, double q, double u, gsl_rng * rand, FILE *fPtr);

double kleinNishinaCrossSection(double energy_ratio);

double scatteredPhotonWeight(double weight, double bias, double optical_depth);
