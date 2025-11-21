void initalizeHotCrossSection(int rank, gsl_rng *rand, FILE *fPtr);

void createHotCrossSection(gsl_rng *rand, FILE *fPtr);

void readHotCrossSection(FILE *fPtr);

int is_dash_line(const char *line);

double calculateTotalThermalCrossSection(double ph_comv, double theta, gsl_rng *rand, FILE *fPtr);

double thermalCrossSectionIntegrand(double x[], size_t dim, void * p);

double boostedCrossSection(double norm_ph_comv, double mu, double gamma);

void display_results (char *title, double result, double error,  FILE *fPtr);

double nonThermalCrossSectionIntegrand(double x[], size_t dim, void * p);

double calculateTotalNonThermalCrossSection(double ph_comv, double gamma_min, double gamma_max, gsl_rng *rand, FILE *fPtr);

void initalizeHotCrossSectionInterp();

double interpolateThermalHotCrossSection(double log_ph_comv_e, double log_theta, FILE *fPtr);

void interpolateSubgroupNonThermalHotCrossSection(double log_ph_comv_e, double *subgroup_interpolated_results, FILE *fPtr);

void cleanupInterpolationData();

void broadcastInterpolationData(int rank);