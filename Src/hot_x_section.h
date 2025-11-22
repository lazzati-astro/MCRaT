// define the extent of the tabulated fluid photon energies normalized by electron rest mass
#define LOG_PH_E_MIN -12.0
#define LOG_PH_E_MAX 6.0
#define N_PH_E 220

// define the extent of the tabulated fluid temperatures normalized by electron rest mass
// for the thermal hot cross section calculation
#define LOG_T_MIN -4.0
#define LOG_T_MAX 4.0
#define N_T 80

//define the thermal modified cross section table filename
#define HOT_THERMAL_X_SECTION_FILE	"thermal_hot_x_section.dat"

#if NONTHERMAL_E_DIST != OFF
    //define the number of lorentz factor intervals that we will calculate the nonthermal hot cross sections for
    #define N_GAMMA 3

    //define the nonthermal modified cross section table filename
    #if NONTHERMAL_E_DIST == POWERLAW
        #define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_pl_hot_x_section.dat"
    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
        #define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_bpl_hot_x_section.dat"
    #else
        #error Unnknown nonthermal electron distribution.
    #endif

#endif


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

#if NONTHERMAL_E_DIST != OFF
void interpolateSubgroupNonThermalHotCrossSection(double log_ph_comv_e, double *subgroup_interpolated_results, FILE *fPtr);
#endif

void cleanupInterpolationData();

void broadcastInterpolationData(int rank);

int checkHotCrossSectionFilesExist(FILE *fPtr);

int validateThermalFile(const char *filename, FILE *fPtr);

#if NONTHERMAL_E_DIST != OFF
int validateNonthermalFile(const char *filename, FILE *fPtr);
#endif