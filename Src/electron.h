//
// Created by Tyler Parsotan on 11/10/25.
//

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr);

double sampleThermalElectron(double temp, gsl_rng * rand, FILE *fPtr);

double samplePowerLaw(double p, double gamma_min, double gamma_max, gsl_rng * rand, FILE *fPtr);

double sampleDoublePowerLaw(double p1, double p2, double gamma_min, double gamma_max, double gamma_break, gsl_rng * rand, FILE *fPtr);
