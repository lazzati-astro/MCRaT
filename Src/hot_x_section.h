void initalizeHotCrossSection(gsl_rng *rand, FILE *fPtr);

double calculateTotalThermalCrossSection(double ph_comv, double theta, gsl_rng *rand, FILE *fPtr);

double thermalCrossSectionIntegrand(double x[], size_t dim, void * p);

double boostedCrossSection(double norm_ph_comv, double mu, double gamma);

void display_results (char *title, double result, double error,  FILE *fPtr);

