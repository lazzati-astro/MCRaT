//
// Created by Tyler Parsotan on 11/18/25.
// This file holds the functions to create/read in the hot cross section lookup table
// This is calculated using the equations outlined in:
// Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009, Astrophysical Journal Supplement, 184, 387 &
// Tomohisa Kawashima et al 2023 ApJ 949 101 &
// Canfield E., Howard W. M. and Liang E. P. 1987 ApJ 323 565
//

#include "mcrat.h"

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

//helper struct to evaluate the double integral
struct double_integral_input_params { double norm_ph_comv; double theta; };

double thermal_table[N_PH_E + 1][N_T + 1];

struct InterpolationData {
    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    double *xa;
    double *ya;
    double *za;
    size_t nx;
    size_t ny;
};

struct InterpolationData global_interp_thermal_data;


#if NONTHERMAL_E_DIST != OFF
    //define the number of lorentz factor intervals that we will calculate the nonthermal hot cross sections for
    #define N_GAMMA 3

    double nonthermal_table[N_PH_E + 1][N_GAMMA];

    //define the nonthermal modified cross section table filename
    #if NONTHERMAL_E_DIST == POWERLAW
        #define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_pl_hot_x_section.dat"
    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
        #define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_bpl_hot_x_section.dat"
    #else
        #error Unnknown nonthermal electron distribution.
    #endif

    struct InterpolationData global_interp_nonthermal_data;


#endif

void initalizeHotCrossSection(int rank, gsl_rng *rand, FILE *fPtr)
{
    int files_exist = 0;

    if (rank == 0)
    {
        // Only rank 0 creates and reads the tables

        // Check if files exist
        files_exist = checkHotCrossSectionFilesExist(fPtr);

        if (!files_exist)
        {
            fprintf(fPtr, "Hot cross section file(s) not found. Creating tables...\n");
            createHotCrossSection(rand, fPtr);
        } else {
            fprintf(fPtr, "Hot cross section file(s) found. Skipping table creation.\n");
        }

        readHotCrossSection(fPtr);
        initalizeHotCrossSectionInterp();
    }

    // Broadcast the interpolation data to all processes
    broadcastInterpolationData(rank);

    fprintf(fPtr, "The hot cross section interpolation has been initialized successfully\n");

    // Test interpolation (all ranks can do this now)
    interpolateThermalHotCrossSection(log10(1e-2), 2.75, fPtr);

    #if NONTHERMAL_E_DIST != OFF
        double test[N_GAMMA];
        interpolateSubgroupNonThermalHotCrossSection(log10(1e-2), test, fPtr);
    #endif
}

void createHotCrossSection(gsl_rng *rand, FILE *fPtr)
{
    int i,j,k;
    double dt=(LOG_T_MAX-LOG_T_MIN)/N_T, dph_e=(LOG_PH_E_MAX-LOG_PH_E_MIN)/N_PH_E;
    double comv_ph_e, theta, gamma_min, gamma_max;
    char xsection_file[STR_BUFFER]="" ;
    FILE *fp;

    for (i = 0; i <= N_PH_E; i++)
    {
        for (j = 0; j <= N_T; j++)
        {
            comv_ph_e = pow(10., LOG_PH_E_MIN + i * dph_e);
            theta = pow(10., LOG_T_MIN + j * dt);
            thermal_table[i][j] = log10(calculateTotalThermalCrossSection(comv_ph_e, theta, rand, fPtr));
            if (isnan(thermal_table[i][j]))
            {
                printf("%d %d %g %g %g\n", i, j, comv_ph_e, theta, thermal_table[i][j]);
                fprintf(fPtr, "%d %d %g %g %g\n", i, j, comv_ph_e, theta, thermal_table[i][j]);
                exit(0);
            }
            else
            {
                printf("%d %d %g %g %g\n", i, j, comv_ph_e, theta, thermal_table[i][j]);
            }
        }
    }

    snprintf(xsection_file,sizeof(xsection_file),"%s%s%s",FILEPATH, MC_PATH,HOT_THERMAL_X_SECTION_FILE);
    fp = fopen(xsection_file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "couldn't write to file\n");
        exit(0);
    }

    //write header to file
    fprintf(fp, "The comoving photon energy and the temperatures are normalized by the electron rest mass\n");
    fprintf(fp, "The calculated hot cross sections are normalized by the thompson cross section.\n");
    fprintf(fp, "Photon index\tTheta Index\tlog10(Comoving Photon Energy)\tlog10(Theta)\tlog10(Hot Cross Section)\n");
    fprintf(fp, "------------------------------------------------\n");


    for (i = 0; i <= N_PH_E; i++)
    {
        for (j = 0; j <= N_T; j++)
        {
            comv_ph_e =  LOG_PH_E_MIN + i * dph_e;
            theta =  LOG_T_MIN + j * dt;
            fprintf(fp, "%d\t%d\t%g\t%g\t%15.10g\n", i, j, comv_ph_e, theta, thermal_table[i][j]);
        }
    }
    fclose(fp);
    fprintf(fPtr, "done.\n\n");

    //set to null incase we are dealing with nonthermal hot cross section after
    fp=NULL;
    memset(&xsection_file[0], 0, sizeof(xsection_file));

    #if NONTHERMAL_E_DIST != OFF
        double dgamma=(log10(GAMMA_MAX)-log10(GAMMA_MIN))/N_GAMMA;

        for (i = 0; i <= N_PH_E; i++)
        {
            //for (j = 0; j <= N_T; j++)
            {
                for (k = 0; k < N_GAMMA; k++)
                {
                    comv_ph_e = pow(10., LOG_PH_E_MIN + i * dph_e);
                    gamma_min = pow(10., log10(GAMMA_MIN) + k * dgamma);
                    gamma_max = pow(10., log10(gamma_min) + dgamma);
                    nonthermal_table[i][k] = log10(calculateTotalNonThermalCrossSection(comv_ph_e, gamma_min, gamma_max,  rand, fPtr));
                    if (isnan(nonthermal_table[i][k]))
                    {
                        fprintf(fPtr, "%d %d %g %g %g %g\n", i, k, comv_ph_e, gamma_min, gamma_max, nonthermal_table[i][k]);
                        printf("%d %d %g %g %g %g\n", i, k, comv_ph_e, gamma_min, gamma_max, nonthermal_table[i][k]);
                        exit(0);
                    }
                    else
                    {
                        printf("%d %d %g %g %g %g\n", i,k, comv_ph_e, gamma_min, gamma_max, nonthermal_table[i][k]);
                    }
                }
            }
        }

        snprintf(xsection_file,sizeof(xsection_file),"%s%s%s",FILEPATH, MC_PATH,HOT_NONTHERMAL_X_SECTION_FILE);

        fp = fopen(xsection_file, "w");
        if (fp == NULL)
        {
            fprintf(stderr, "couldn't write to file\n");
            exit(0);
        }

        #if NONTHERMAL_E_DIST == POWERLAW
            fprintf(fp, "This file is produced for a Powerlaw electron distribution with:\n");
            fprintf(fp, "gamma_min: %g\ngamma_max: %g\npowerlaw index %g\n", GAMMA_MIN, GAMMA_MAX, POWERLAW_INDEX);
        #else
            fprintf(fp, "This file is produced for a Broken Powerlaw electron distribution with:\n");
            fprintf(fp, "gamma_min: %g\ngamma_max: %g\ngamma_break: %g\n powerlaw index 1 %g\npowerlaw index 2 %g\n", GAMMA_MIN, GAMMA_MAX, GAMMA_BREAK, POWERLAW_INDEX_1, POWERLAW_INDEX_2);
        #endif
        fprintf(fp, "The comoving photon energy is normalized by the electron rest mass\n");
        fprintf(fp, "The calculated hot cross sections are normalized by the thompson cross section.\n");
        fprintf(fp, "Photon index\telectron gamma Index\tlog10(Comoving Photon Energy)\tlog10(electron_gamma_min)\tlog10(electron_gamma_max)\tlog10(Hot Cross Section)\n");
        fprintf(fp, "------------------------------------------------\n");

        for (i = 0; i <= N_PH_E; i++)
        {
            //for (j = 0; j <= N_T; j++)
            {
                for (k = 0; k < N_GAMMA; k++)
                {
                    comv_ph_e =  LOG_PH_E_MIN + i * dph_e;
                    gamma_min = log10(GAMMA_MIN) + k * dgamma;
                    gamma_max = gamma_min + dgamma;
                    fprintf(fp, "%d\t%d\t%g\t%g\t%g\t%15.10g\n", i, k, comv_ph_e, gamma_min, gamma_max, nonthermal_table[i][k]);
                }
            }
        }
        fclose(fp);
        fprintf(stderr, "done.\n\n");

    #endif
}

void readHotCrossSection(FILE *fPtr)
{
    int i, j, read_error, parsed;
    double comv_ph_e, theta, gamma_min, gamma_max, value;
    char xsection_file[STR_BUFFER]="", line[STR_BUFFER] ;
    FILE *fp;

    snprintf(xsection_file,sizeof(xsection_file),"%s%s%s",FILEPATH, MC_PATH,HOT_THERMAL_X_SECTION_FILE);

    fprintf(fPtr, "Reading thermal hot cross section data from %s...\n", xsection_file);

    //first read in the thermal electron hot cross section
    fp = fopen(xsection_file, "r");

    if (fp == NULL)
    {
        fprintf(fPtr, "Couldn't read file\n");
        exit(0);
    }

    // Skip header lines until line with only dashes is found
    fgets(line, sizeof(line), fp);
    line[strcspn(line, "\r\n")] = 0;
    while (fgets(line, sizeof(line), fp) && !is_dash_line(line))
    {
        // Remove trailing newline
        line[strcspn(line, "\r\n")] = 0;
    }

    // Read data rows
    while (fgets(line, sizeof(line), fp))
    {
        parsed = sscanf(line, "%d\t%d\t%lf\t%lf\t%lf", &i, &j, &comv_ph_e, &theta, &value);
        if (i >=0  && i <= N_PH_E && j >= 0 && j <= N_T)
        {
            thermal_table[i][j] = value;
        }
        else
        {
            fprintf(fPtr, "The bounds of the input file exceed what MCRaT has been compiled with.\n");
            exit(0);
        }
    }

    fclose(fp);
    //set to null incase we are dealing with nonthermal hot cross section after
    fp=NULL;
    memset(&xsection_file[0], 0, sizeof(xsection_file));


    #if NONTHERMAL_E_DIST != OFF
        //TODO: add check for whether the nonthermal electron distribution valuess match in teh header match the compiler macro values

        snprintf(xsection_file,sizeof(xsection_file),"%s%s%s",FILEPATH, MC_PATH,HOT_NONTHERMAL_X_SECTION_FILE);

        fprintf(fPtr, "Reading nonthermal hot cross section data from %s...\n", xsection_file);

        //first read in the thermal electron hot cross section
        fp = fopen(xsection_file, "r");

        if (fp == NULL)
        {
            fprintf(fPtr, "Couldn't read file\n");
            exit(0);
        }

        // Skip header lines until line with only dashes is found
        fgets(line, sizeof(line), fp);
        line[strcspn(line, "\r\n")] = 0;
        while (fgets(line, sizeof(line), fp) && !is_dash_line(line))
        {
            // Remove trailing newline
            line[strcspn(line, "\r\n")] = 0;
        }

        // Read data rows
        while (fgets(line, sizeof(line), fp))
        {
            parsed = sscanf(line, "%d\t%d\t%lf\t%lf\t%lf\t%lf", &i, &j, &comv_ph_e, &gamma_min, &gamma_max, &value);
            if (i >=0  && i <= N_PH_E && j >= 0 && j < N_GAMMA)
            {
                nonthermal_table[i][j] = value;
            }
            else
            {
                fprintf(fPtr, "The bounds of the non-thermal input file exceed what MCRaT has been compiled with.\n");
                exit(0);
            }
        }

        fclose(fp);

    #endif
}

int is_dash_line(const char *line)
{
    int found_dash = 0;
    for (int i = 0; line[i]; i++)
    {
        if (!isspace((unsigned char)line[i]))
        {
            if (line[i] != '-')
            {
                return 0;
            }
            found_dash = 1;
        }
    }
    return found_dash;
}

double calculateTotalThermalCrossSection(double ph_comv, double theta, gsl_rng *rand, FILE *fPtr)
{
    /*
    Here we perform a double integral over gamma and mu, given the photon energy (normalized by the electron rest mass)
    in the fluid rest frame and the fluid temperature. This integral is done using the maxwell juttner distribution of
    thermal electrons. The output is normalized by the thompson cross section and includes the 0.5 factor that is
    included in the full definition of the hot cross section.
    */
    double result=0, error=0;
    double xl[2] = { 1, -1 };
    double xu[2] = { 1. + 12 * theta, 1 };
    struct double_integral_input_params params = {ph_comv, theta };

    //
    if (theta < pow(10, LOG_T_MIN) && ph_comv < pow(10, LOG_PH_E_MIN))
        return 1;
    if (theta < pow(10, LOG_T_MIN))
        return (kleinNishinaCrossSection(ph_comv));

    gsl_monte_function F;
    F.f = &thermalCrossSectionIntegrand;
    F.dim = 2;
    F.params = &params;

    size_t calls = 500000;

    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&F, xl, xu, 2, calls, rand, s, &result, &error);
    gsl_monte_plain_free (s);

    display_results ("plain", result, error, fPtr);

    return 0.5*result;
}

double thermalCrossSectionIntegrand(double x[], size_t dim, void * p)
{
    double result=0;
    struct double_integral_input_params * fp = (struct double_integral_input_params *)p;
    double gamma=x[0], mu=x[1];

    result = singleMaxwellJuttner(gamma, fp->theta)*boostedCrossSection(fp->norm_ph_comv, mu, gamma);

    return result;
}

double boostedCrossSection(double norm_ph_comv, double mu, double gamma)
{
    /*
        Calculates the KN cross section normalized by the thompson cross section in the electron rest frame
        This takes the photon's comobetaing energy (in the fluid frame) normalized by the electron rest mass, norm_ph_comv,
        the cosine of the angle between the photon and the electron in the electron rest frame, mu,
        and the lorentz factor of the electron in the fluid frame, gamma
    */
    double norm_ph_e=0, result=0, beta=0;

    /* energy in electron rest frame */
    beta = sqrt(gamma * gamma - 1.) / gamma;
    norm_ph_e = norm_ph_comv * gamma * (1. - mu * beta);

    result = kleinNishinaCrossSection(norm_ph_e) * (1. - mu * beta);

    if (result > 2)
    {
        fprintf(stderr, "norm_ph_comv,mu,gamma: %g %g %g\n", norm_ph_comv, mu, gamma);
        fprintf(stderr, "beta,norm_ph_e, result: %g %g %g\n", beta, norm_ph_e, result);
        fprintf(stderr, "kn: %g %g %g\n", beta, norm_ph_e, result);
    }

    if (isnan(result))
    {
        fprintf(stderr, "isnan: %g %g %g\n", norm_ph_comv, mu, gamma);
        exit(0);
    }

    return result;
}

void display_results (char *title, double result, double error,  FILE *fPtr)
{
    fprintf (fPtr, "%s ==================\n", title);
    fprintf (fPtr, "result = % .6f\n", result);
    fprintf (fPtr, "sigma  = % .6f\n", error);
    //printf ("exact  = % .6f\n", exact);
    //printf ("error  = % .6f = %.2g sigma\n", result - exact, fabs (result - exact) / error);
}

double nonThermalCrossSectionIntegrand(double x[], size_t dim, void * p)
{
    double result=0;
    struct double_integral_input_params * fp = (struct double_integral_input_params *)p;
    double gamma=x[0], mu=x[1];
    #if NONTHERMAL_E_DIST != OFF

        #if NONTHERMAL_E_DIST == POWERLAW
            result = singleElectronPowerLaw(gamma, POWERLAW_INDEX, GAMMA_MIN, GAMMA_MAX);
        #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
            result = singleElectronBrokenPowerLaw(gamma, POWERLAW_INDEX_1, POWERLAW_INDEX_2, GAMMA_MIN, GAMMA_MAX, GAMMA_BREAK);
        #else
            #error Unnknown nonthermal electron distribution.
        #endif
    #endif

    result*=boostedCrossSection(fp->norm_ph_comv, mu, gamma);

    return result;
}

double calculateTotalNonThermalCrossSection(double ph_comv, double gamma_min, double gamma_max, gsl_rng *rand, FILE *fPtr)
{
    /*
    Here we perform a double integral over gamma and mu, given the photon energy (normalized by the electron rest mass)
    in the fluid rest frame and the fluid temperature. This integral is done using the maxwell juttner distribution of
    thermal electrons. The output is normalized by the thompson cross section and includes the 0.5 factor that is
    included in the full definition of the hot cross section.
    */
    double result=0, error=0;
    double xl[2] = { gamma_min, -1 };
    double xu[2] = { gamma_max, 1 };
    struct double_integral_input_params params = {ph_comv, 0.0 };

    gsl_monte_function F;
    F.f = &nonThermalCrossSectionIntegrand;
    F.dim = 2;
    F.params = &params;

    size_t calls = 1000000;

    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&F, xl, xu, 2, calls, rand, s, &result, &error);
    gsl_monte_plain_free (s);

    display_results ("plain", result, error, fPtr);

    return 0.5*result;
}

void initalizeHotCrossSectionInterp()
{
    int i=0, j=0;
    double dt=(LOG_T_MAX-LOG_T_MIN)/N_T, dph_e=(LOG_PH_E_MAX-LOG_PH_E_MIN)/N_PH_E, dgamma=0;
    // Dynamically allocate arrays for the interpolation struct
    double *comv_ph_grid = malloc((N_PH_E + 1) * sizeof(double));
    double *theta_grid = malloc((N_T + 1) * sizeof(double));
    double *data_grid = malloc((N_PH_E + 1) * (N_T + 1) * sizeof(double));

    for (i = 0; i <= N_PH_E; i++)
    {
        comv_ph_grid[i] = LOG_PH_E_MIN + i * dph_e;
        printf("comv_ph_grid[%d] = %g\n", i, comv_ph_grid[i]);
    }

    for (i = 0; i <= N_T; i++)
    {
        theta_grid[i] = LOG_T_MIN + i * dt;
        printf("theta_grid[%d] = %g\n", i, theta_grid[i]);
    }

    for (i = 0; i <= N_PH_E; i++)
    {
        for (j = 0; j <= N_T; j++)
        {
            data_grid[j*(N_PH_E+1)+i] = thermal_table[i][j];
        }
    }

    global_interp_thermal_data.nx = N_PH_E+1;
    global_interp_thermal_data.ny = N_T+1;
    global_interp_thermal_data.xa = comv_ph_grid;
    global_interp_thermal_data.ya = theta_grid;
    global_interp_thermal_data.za = data_grid;

    // The gsl_spline2d high-level interface stores the data arrays internally
    global_interp_thermal_data.spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, N_PH_E+1, N_T+1);
    global_interp_thermal_data.xacc = gsl_interp_accel_alloc();
    global_interp_thermal_data.yacc = gsl_interp_accel_alloc();

    // Initialize the spline with the data
    gsl_spline2d_init(global_interp_thermal_data.spline, global_interp_thermal_data.xa, global_interp_thermal_data.ya, global_interp_thermal_data.za, global_interp_thermal_data.nx, global_interp_thermal_data.ny);

    #if NONTHERMAL_E_DIST != OFF
        dgamma=(log10(GAMMA_MAX)-log10(GAMMA_MIN))/N_GAMMA;
        double gamma_min=0, gamma_max=0;
        double *gamma_grid = malloc(N_GAMMA * sizeof(double));
        double *nonthermal_data_grid = malloc(N_GAMMA * (N_PH_E + 1) * sizeof(double));

        for (i = 0; i < N_GAMMA; i++)
        {
            gamma_min = log10(GAMMA_MIN) + i * dgamma;
            gamma_max = gamma_min + dgamma;
            gamma_grid[i] = 0.5*(gamma_min+gamma_max);
            printf("gamma_grid[%d] = %g\n", i, gamma_grid[i]);
        }


        for (i = 0; i <= N_PH_E; i++)
        {
            for (j = 0; j < N_GAMMA; j++)
            {
                nonthermal_data_grid[j*(N_PH_E+1)+i] = nonthermal_table[i][j];
            }
        }

        global_interp_nonthermal_data.nx = N_PH_E+1;
        global_interp_nonthermal_data.ny = N_GAMMA;
        global_interp_nonthermal_data.xa = comv_ph_grid;
        global_interp_nonthermal_data.ya = gamma_grid;
        global_interp_nonthermal_data.za = nonthermal_data_grid;

        // The gsl_spline2d high-level interface stores the data arrays internally
        global_interp_nonthermal_data.spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, N_PH_E+1, N_GAMMA);
        global_interp_nonthermal_data.xacc = gsl_interp_accel_alloc();
        global_interp_nonthermal_data.yacc = gsl_interp_accel_alloc();

        // Initialize the spline with the data
        gsl_spline2d_init(global_interp_nonthermal_data.spline, global_interp_nonthermal_data.xa, global_interp_nonthermal_data.ya, global_interp_nonthermal_data.za, global_interp_nonthermal_data.nx, global_interp_nonthermal_data.ny);
    #endif

}

//interpolation checked with python interpolation of the same hot coss section table
double interpolateThermalHotCrossSection(double log_ph_comv_e, double log_theta, FILE *fPtr)
{
    double result=0;
    // Access global_interp_data fields
    result=gsl_spline2d_eval(global_interp_thermal_data.spline, log_ph_comv_e, log_theta, global_interp_thermal_data.xacc, global_interp_thermal_data.yacc);
    fprintf(fPtr, "Thermal: %g %g %g\n", log_ph_comv_e, log_theta, result);
    return result;
}

void interpolateSubgroupNonThermalHotCrossSection(double log_ph_comv_e, double *subgroup_interpolated_results, FILE *fPtr)
{
    // iterate over the subgroups to get the nonthermal cross sections and save them to the pointer array
    int i=0;
    double results[N_GAMMA];

    for (i=0;i<global_interp_nonthermal_data.ny;i++)
    {
        results[i]=gsl_spline2d_eval(global_interp_nonthermal_data.spline, log_ph_comv_e, global_interp_nonthermal_data.ya[i], global_interp_nonthermal_data.xacc, global_interp_nonthermal_data.yacc);
        fprintf(fPtr, "Non-thermal: %g %g %g\n", log_ph_comv_e, global_interp_nonthermal_data.ya[i], results[i]);
    }

    //todo: make sure that the pointer has enough space allocated
    for (i=0;i<global_interp_nonthermal_data.ny;i++)
    {
        *(subgroup_interpolated_results+i)=results[i];
    }
}

void cleanupInterpolationData()
{
    gsl_spline2d_free(global_interp_thermal_data.spline);
    gsl_interp_accel_free(global_interp_thermal_data.xacc);
    gsl_interp_accel_free(global_interp_thermal_data.yacc);

    // Free dynamically allocated arrays
    free(global_interp_thermal_data.xa);
    free(global_interp_thermal_data.ya);
    free(global_interp_thermal_data.za);

    #if NONTHERMAL_E_DIST != OFF
        gsl_spline2d_free(global_interp_nonthermal_data.spline);
        gsl_interp_accel_free(global_interp_nonthermal_data.xacc);
        gsl_interp_accel_free(global_interp_nonthermal_data.yacc);

        // Free nonthermal arrays (xa is shared with thermal, so don't free twice)
        free(global_interp_nonthermal_data.ya);
        free(global_interp_nonthermal_data.za);
    #endif
}

void broadcastInterpolationData(int rank)
{
    int i, j;
    double dt = (LOG_T_MAX - LOG_T_MIN) / N_T;
    double dph_e = (LOG_PH_E_MAX - LOG_PH_E_MIN) / N_PH_E;

    // Broadcast thermal table data
    MPI_Bcast(thermal_table, (N_PH_E + 1) * (N_T + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        // Non-root processes need to allocate and initialize their interpolation structures

        // Allocate arrays for grid points and data
        double *comv_ph_grid = malloc((N_PH_E + 1) * sizeof(double));
        double *theta_grid = malloc((N_T + 1) * sizeof(double));
        double *data_grid = malloc((N_PH_E + 1) * (N_T + 1) * sizeof(double));

        // Reconstruct grid points
        for (i = 0; i <= N_PH_E; i++)
        {
            comv_ph_grid[i] = LOG_PH_E_MIN + i * dph_e;
        }

        for (i = 0; i <= N_T; i++)
        {
            theta_grid[i] = LOG_T_MIN + i * dt;
        }

        // Copy thermal table data to data_grid
        for (i = 0; i <= N_PH_E; i++)
        {
            for (j = 0; j <= N_T; j++)
            {
                data_grid[j * (N_PH_E + 1) + i] = thermal_table[i][j];
            }
        }

        // Initialize interpolation structure
        global_interp_thermal_data.nx = N_PH_E + 1;
        global_interp_thermal_data.ny = N_T + 1;
        global_interp_thermal_data.xa = comv_ph_grid;
        global_interp_thermal_data.ya = theta_grid;
        global_interp_thermal_data.za = data_grid;

        global_interp_thermal_data.spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, N_PH_E + 1, N_T + 1);
        global_interp_thermal_data.xacc = gsl_interp_accel_alloc();
        global_interp_thermal_data.yacc = gsl_interp_accel_alloc();

        gsl_spline2d_init(global_interp_thermal_data.spline,
                          global_interp_thermal_data.xa,
                          global_interp_thermal_data.ya,
                          global_interp_thermal_data.za,
                          global_interp_thermal_data.nx,
                          global_interp_thermal_data.ny);
    }

    #if NONTHERMAL_E_DIST != OFF
        // Broadcast nonthermal table data
        MPI_Bcast(nonthermal_table, (N_PH_E + 1) * N_GAMMA, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank != 0)
        {
            double dgamma = (log10(GAMMA_MAX) - log10(GAMMA_MIN)) / N_GAMMA;
            double gamma_min, gamma_max;

            // Allocate arrays
            double *comv_ph_grid_nt = malloc((N_PH_E + 1) * sizeof(double));
            double *gamma_grid = malloc(N_GAMMA * sizeof(double));
            double *nonthermal_data_grid = malloc((N_PH_E + 1) * N_GAMMA * sizeof(double));

            // Reconstruct photon energy grid
            for (i = 0; i <= N_PH_E; i++)
            {
                comv_ph_grid_nt[i] = LOG_PH_E_MIN + i * dph_e;
            }

            // Reconstruct gamma grid
            for (i = 0; i < N_GAMMA; i++)
            {
                gamma_min = log10(GAMMA_MIN) + i * dgamma;
                gamma_max = gamma_min + dgamma;
                gamma_grid[i] = 0.5 * (gamma_min + gamma_max);
            }

            // Copy nonthermal table data
            for (i = 0; i <= N_PH_E; i++)
            {
                for (j = 0; j < N_GAMMA; j++)
                {
                    nonthermal_data_grid[j * (N_PH_E + 1) + i] = nonthermal_table[i][j];
                }
            }

            // Initialize nonthermal interpolation structure
            global_interp_nonthermal_data.nx = N_PH_E + 1;
            global_interp_nonthermal_data.ny = N_GAMMA;
            global_interp_nonthermal_data.xa = comv_ph_grid_nt;
            global_interp_nonthermal_data.ya = gamma_grid;
            global_interp_nonthermal_data.za = nonthermal_data_grid;

            global_interp_nonthermal_data.spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, N_PH_E + 1, N_GAMMA);
            global_interp_nonthermal_data.xacc = gsl_interp_accel_alloc();
            global_interp_nonthermal_data.yacc = gsl_interp_accel_alloc();

            gsl_spline2d_init(global_interp_nonthermal_data.spline,
                              global_interp_nonthermal_data.xa,
                              global_interp_nonthermal_data.ya,
                              global_interp_nonthermal_data.za,
                              global_interp_nonthermal_data.nx,
                              global_interp_nonthermal_data.ny);
        }
    #endif

    // Ensure all processes are synchronized
    MPI_Barrier(MPI_COMM_WORLD);
}

int checkHotCrossSectionFilesExist(FILE *fPtr)
{
    char thermal_file[STR_BUFFER] = "";
    int thermal_valid = 0;
    int nonthermal_valid = 1; // Default to true if not needed

    // Check and validate thermal file
    snprintf(thermal_file, sizeof(thermal_file), "%s%s%s",
             FILEPATH, MC_PATH, HOT_THERMAL_X_SECTION_FILE);

    thermal_valid = validateThermalFile(thermal_file, fPtr);

    #if NONTHERMAL_E_DIST != OFF
        char nonthermal_file[STR_BUFFER] = "";

        snprintf(nonthermal_file, sizeof(nonthermal_file), "%s%s%s",
                 FILEPATH, MC_PATH, HOT_NONTHERMAL_X_SECTION_FILE);

        nonthermal_valid = validateNonthermalFile(nonthermal_file, fPtr);
    #endif

    return thermal_valid && nonthermal_valid;
}

int validateThermalFile(const char *filename, FILE *fPtr)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(fPtr, "Thermal file does not exist: %s\n", filename);
        return 0;
    }

    char line[STR_BUFFER];
    int i, j;
    double comv_ph_e, theta, value;
    int data_lines = 0;
    int in_data_section = 0;
    int grid_valid = 1;

    // Calculate expected grid values
    double dt = (LOG_T_MAX - LOG_T_MIN) / N_T;
    double dph_e = (LOG_PH_E_MAX - LOG_PH_E_MIN) / N_PH_E;
    double expected_comv_ph, expected_theta;
    double tolerance = 1e-9; // Tolerance for floating point comparison of the tabulated values
    double tolerance_grid=1e-3;    //tolerance for checking the grid values that were used to create the table

    // Skip header lines until dashes
    while (fgets(line, sizeof(line), fp))
    {
        line[strcspn(line, "\r\n")] = 0;
        if (is_dash_line(line))
        {
            in_data_section = 1;
            break;
        }
    }

    if (!in_data_section)
    {
        fprintf(fPtr, "Thermal file invalid: header format incorrect in %s\n", filename);
        fclose(fp);
        return 0;
    }

    // Read and validate data lines
    while (fgets(line, sizeof(line), fp))
        {
        int parsed = sscanf(line, "%d\t%d\t%lf\t%lf\t%lf", &i, &j, &comv_ph_e, &theta, &value);

        if (parsed != 5)
        {
            continue; // Skip malformed lines
        }

        // Check indices are in valid range
        if (i < 0 || i > N_PH_E || j < 0 || j > N_T)
        {
            fprintf(fPtr, "Thermal file invalid: indices out of range [%d, %d] in %s\n",
                    i, j, filename);
            grid_valid = 0;
            break;
        }

        // Calculate expected grid values
        expected_comv_ph = LOG_PH_E_MIN + i * dph_e;
        expected_theta = LOG_T_MIN + j * dt;

        // Validate grid values
        if (fabs(comv_ph_e - expected_comv_ph) > tolerance_grid)
        {
            fprintf(fPtr, "Thermal file invalid: comv_ph_e mismatch at [%d, %d]. "
                    "Expected %g, found %g with Delta: %e in %s\n",
                    i, j, expected_comv_ph, comv_ph_e, fabs(comv_ph_e - expected_comv_ph), filename);
            grid_valid = 0;
            break;
        }

        if (fabs(theta - expected_theta) > tolerance_grid)
        {
            fprintf(fPtr, "Thermal file invalid: theta mismatch at [%d, %d]. "
                    "Expected %g, found %g with Delta: %e in %s\n",
                    i, j, expected_theta, theta, fabs(theta - expected_theta), filename);
            grid_valid = 0;
            break;
        }

        data_lines++;
    }

    fclose(fp);

    // Expected number of data lines
    int expected_lines = (N_PH_E + 1) * (N_T + 1);

    if (data_lines != expected_lines)
    {
        fprintf(fPtr, "Thermal file invalid: %s (expected %d lines, found %d)\n",
                filename, expected_lines, data_lines);
        return 0;
    }

    if (!grid_valid)
    {
        return 0;
    }

    fprintf(fPtr, "Thermal file validated successfully: %s (%d lines)\n", filename, data_lines);
    return 1;
}

#if NONTHERMAL_E_DIST != OFF
    int validateNonthermalFile(const char *filename, FILE *fPtr)
    {
        FILE *fp = fopen(filename, "r");
        if (fp == NULL)
        {
            fprintf(fPtr, "Nonthermal file does not exist: %s\n", filename);
            return 0;
        }

        char line[STR_BUFFER];
        int i, j;
        double comv_ph_e, gamma_min, gamma_max, value;
        int data_lines = 0;
        int in_data_section = 0;
        int header_valid = 0;
        int grid_valid = 1;

        // Calculate expected grid values
        double dph_e = (LOG_PH_E_MAX - LOG_PH_E_MIN) / N_PH_E;
        double dgamma = (log10(GAMMA_MAX) - log10(GAMMA_MIN)) / N_GAMMA;
        double expected_comv_ph, expected_gamma_min, expected_gamma_max;
        double tolerance = 1e-9; // Tolerance for floating point comparison of the tabulated values
        double tolerance_grid=1e-3;    //tolerance for checking the grid values that were used to create the table


        // Read header and validate parameters
        while (fgets(line, sizeof(line), fp))
        {
            line[strcspn(line, "\r\n")] = 0;

            if (is_dash_line(line))
            {
                in_data_section = 1;
                break;
            }

            // Validate header parameters match compilation settings
            #if NONTHERMAL_E_DIST == POWERLAW
                double file_gamma_min, file_gamma_max, file_index;

                if (sscanf(line, "gamma_min: %lf", &file_gamma_min) == 1)
                {
                    if (fabs(file_gamma_min - GAMMA_MIN) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: GAMMA_MIN mismatch. "
                                "Expected %g, found %g in %s\n",
                                GAMMA_MIN, file_gamma_min, filename);
                        fclose(fp);
                        return 0;
                    }
                }

                if (sscanf(line, "gamma_max: %lf", &file_gamma_max) == 1)
                {
                    if (fabs(file_gamma_max - GAMMA_MAX) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: GAMMA_MAX mismatch. "
                                "Expected %g, found %g in %s\n",
                                GAMMA_MAX, file_gamma_max, filename);
                        fclose(fp);
                        return 0;
                    }
                }

                if (sscanf(line, "powerlaw index %lf", &file_index) == 1)
                {
                    if (fabs(file_index - POWERLAW_INDEX) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: POWERLAW_INDEX mismatch. "
                                "Expected %g, found %g in %s\n",
                                POWERLAW_INDEX, file_index, filename);
                        fclose(fp);
                        return 0;
                    }
                    header_valid = 1;
                }

            #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
                double file_gamma_min, file_gamma_max, file_gamma_break;
                double file_index1, file_index2;

                if (sscanf(line, "gamma_min: %lf", &file_gamma_min) == 1)
                {
                    if (fabs(file_gamma_min - GAMMA_MIN) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: GAMMA_MIN mismatch. "
                                "Expected %g, found %g in %s\n",
                                GAMMA_MIN, file_gamma_min, filename);
                        fclose(fp);
                        return 0;
                    }
                }

                if (sscanf(line, "gamma_max: %lf", &file_gamma_max) == 1)
                {
                    if (fabs(file_gamma_max - GAMMA_MAX) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: GAMMA_MAX mismatch. "
                                "Expected %g, found %g in %s\n",
                                GAMMA_MAX, file_gamma_max, filename);
                        fclose(fp);
                        return 0;
                    }
                }

                if (sscanf(line, "gamma_break: %lf", &file_gamma_break) == 1)
                {
                    if (fabs(file_gamma_break - GAMMA_BREAK) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: GAMMA_BREAK mismatch. "
                                "Expected %g, found %g in %s\n",
                                GAMMA_BREAK, file_gamma_break, filename);
                        fclose(fp);
                        return 0;
                    }
                }

                if (sscanf(line, "powerlaw index 1 %lf", &file_index1) == 1)
                {
                    if (fabs(file_index1 - POWERLAW_INDEX_1) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: POWERLAW_INDEX_1 mismatch. "
                                "Expected %g, found %g in %s\n",
                                POWERLAW_INDEX_1, file_index1, filename);
                        fclose(fp);
                        return 0;
                    }
                }

                if (sscanf(line, "powerlaw index 2 %lf", &file_index2) == 1)
                {
                    if (fabs(file_index2 - POWERLAW_INDEX_2) > tolerance)
                    {
                        fprintf(fPtr, "Nonthermal file invalid: POWERLAW_INDEX_2 mismatch. "
                                "Expected %g, found %g in %s\n",
                                POWERLAW_INDEX_2, file_index2, filename);
                        fclose(fp);
                        return 0;
                    }
                    header_valid = 1;
                }
            #endif
        }

        if (!in_data_section)
        {
            fprintf(fPtr, "Nonthermal file invalid: header format incorrect in %s\n", filename);
            fclose(fp);
            return 0;
        }

        if (!header_valid)
        {
            fprintf(fPtr, "Nonthermal file invalid: required parameters not found in header of %s\n",
                    filename);
            fclose(fp);
            return 0;
        }

        // Read and validate data lines
        while (fgets(line, sizeof(line), fp))
        {
            int parsed = sscanf(line, "%d\t%d\t%lf\t%lf\t%lf\t%lf",
                               &i, &j, &comv_ph_e, &gamma_min, &gamma_max, &value);

            if (parsed != 6) {
                continue; // Skip malformed lines
            }

            // Check indices are in valid range
            if (i < 0 || i > N_PH_E || j < 0 || j >= N_GAMMA)
            {
                fprintf(fPtr, "Nonthermal file invalid: indices out of range [%d, %d] in %s\n",
                        i, j, filename);
                grid_valid = 0;
                break;
            }

            // Calculate expected grid values
            expected_comv_ph = LOG_PH_E_MIN + i * dph_e;
            expected_gamma_min = log10(GAMMA_MIN) + j * dgamma;
            expected_gamma_max = expected_gamma_min + dgamma;

            // Validate comv_ph_e grid
            if (fabs(comv_ph_e - expected_comv_ph) > tolerance_grid)
            {
                fprintf(fPtr, "Nonthermal file invalid: comv_ph_e mismatch at [%d, %d]. "
                        "Expected %g, found %g with Delta: %e in %s\n",
                        i, j, expected_comv_ph, comv_ph_e, fabs(comv_ph_e - expected_comv_ph), filename);
                grid_valid = 0;
                break;
            }

            // Validate gamma_min grid
            if (fabs(gamma_min - expected_gamma_min) > tolerance_grid)
            {
                fprintf(fPtr, "Nonthermal file invalid: gamma_min mismatch at [%d, %d]. "
                        "Expected %g, found %g with Delta: %e, in %s\n",
                        i, j, expected_gamma_min, gamma_min,fabs(gamma_min - expected_gamma_min), filename);
                grid_valid = 0;
                break;
            }

            // Validate gamma_max grid
            if (fabs(gamma_max - expected_gamma_max) > tolerance_grid)
            {
                fprintf(fPtr, "Nonthermal file invalid: gamma_max mismatch at [%d, %d]. "
                        "Expected %g, found %g with Delta: %e in %s\n",
                        i, j, expected_gamma_max, gamma_max, fabs(gamma_max - expected_gamma_max), filename);
                grid_valid = 0;
                break;
            }

            data_lines++;
        }

        fclose(fp);

        // Expected number of data lines
        int expected_lines = (N_PH_E + 1) * N_GAMMA;

        if (data_lines != expected_lines)
        {
            fprintf(fPtr, "Nonthermal file invalid: %s (expected %d lines, found %d)\n",
                    filename, expected_lines, data_lines);
            return 0;
        }

        if (!grid_valid)
        {
            return 0;
        }

        fprintf(fPtr, "Nonthermal file validated successfully: %s (%d lines)\n",
                filename, data_lines);
        return 1;
    }
#endif



