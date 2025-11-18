//
// Created by Tyler Parsotan on 11/10/25.
//

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr);

double sampleElectronTheta(double beta, gsl_rng * rand, FILE *fPtr);

double sampleThermalElectron(double temp, gsl_rng * rand, FILE *fPtr);

double sampleNonthermalElectron(double p, gsl_rng * rand, FILE *fPtr);

double samplePowerLaw(double p, double gamma_min, double gamma_max, gsl_rng * rand, FILE *fPtr);

double sampleBrokenPowerLaw(double p1, double p2, double gamma_min, double gamma_max, double gamma_break, gsl_rng * rand, FILE *fPtr);

double brokenPowerLawNorm(double p1, double p2, double gamma_min, double gamma_max, double gamma_break);

double singleElectronBrokenPowerLaw(double x, double p1, double p2, double gamma_min, double gamma_max, double gamma_break);

void arrayElectronBrokenPowerLaw(const double *x, double *y, int n_points, double p1, double p2, double gamma_min, double gamma_max, double gamma_break);

double powerLawNorm(double p, double gamma_min, double gamma_max);

double singleElectronPowerLaw(double x, double p, double gamma_min, double gamma_max);

void arrayElectronPowerLaw(const double *x, double *y, int n_points, double p, double gamma_min, double gamma_max);

double singleMaxwellJuttner(double gamma, double theta);