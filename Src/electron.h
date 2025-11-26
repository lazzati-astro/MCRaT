//
// Created by Tyler Parsotan on 11/10/25.
//

int generateSingleElectron(double *el_p, double temp, double *ph_p, struct photon *ph, gsl_rng * rand, FILE *fPtr);

void singleThermalElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr);

void singleNonThermalElectron(double *el_p, double *ph_p, double gamma_min, double gamma_max, gsl_rng * rand, FILE *fPtr);

void rotateElectron(double *el_p, double *ph_p, FILE *fPtr);

double sampleElectronTheta(double beta, gsl_rng * rand, FILE *fPtr);

double sampleThermalElectron(double temp, gsl_rng * rand, FILE *fPtr);

double sampleNonthermalElectron(gsl_rng * rand, FILE *fPtr);

double samplePowerLaw(double p, double gamma_min, double gamma_max, gsl_rng * rand, FILE *fPtr);

double sampleBrokenPowerLaw(double p1, double p2, double gamma_min, double gamma_max, double gamma_break, gsl_rng * rand, FILE *fPtr);

double brokenPowerLawNorm(double p1, double p2, double gamma_min, double gamma_max, double gamma_break);

double singleElectronBrokenPowerLaw(double x, double p1, double p2, double gamma_min, double gamma_max, double gamma_break);

void arrayElectronBrokenPowerLaw(const double *x, double *y, int n_points, double p1, double p2, double gamma_min, double gamma_max, double gamma_break);

double powerLawNorm(double p, double gamma_min, double gamma_max);

double singleElectronPowerLaw(double x, double p, double gamma_min, double gamma_max);

void arrayElectronPowerLaw(const double *x, double *y, int n_points, double p, double gamma_min, double gamma_max);

double singleMaxwellJuttner(double gamma, double theta);

double nonThermalElectronDistIntegrand(double x, void * params);

double calculateNormPowerLawEnergyDens(double p, double gamma_min, double gamma_max);

double calculateNormBrokenPowerLawEnergyDens(double p1, double p2, double gamma_min, double gamma_max, double gamma_break);


#if NONTHERMAL_E_DIST != OFF
    void calculateElectronDistSubgroupDens(double *subgroup_dens, FILE *fPtr);

    void calculateNonthermalElectronDens(struct hydro_dataframe *hydro_data, FILE *fPtr);
#endif