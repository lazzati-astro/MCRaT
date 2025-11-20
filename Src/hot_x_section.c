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

//define the number of lorentz factor intervals that we will calculate the nonthermal hot cross sections for
#define N_GAMMA 3

//helper struct to evaluate the double integral
struct double_integral_input_params { double norm_ph_comv; double theta };

double thermal_table[N_PH_E + 1][N_T + 1];

#ifdef NONTHERMAL_E_DIST
    double nonthermal_table[N_PH_E + 1][N_GAMMA];

    //define the nonthermal modified cross section table filename
    #if NONTHERMAL_E_DIST == POWERLAW
        #define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_pl_hot_x_section.dat"
    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
        #define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_bpl_hot_x_section.dat"
    #else
        #error Unnknown nonthermal electron distribution.
    #endif

#endif

void initalizeHotCrossSection(gsl_rng *rand, FILE *fPtr)
{
    createHotCrossSection(rand, fPtr)

    readHotCrossSection(fPtr)
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

    #ifdef NONTHERMAL_E_DIST
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
                    fprintf(fp, "%d %d %g %g %g %15.10g\n", i, k, comv_ph_e, gamma_min, gamma_max, nonthermal_table[i][k]);
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
        if (i >=0  && i < N_PH_E && j >= 0 && j < N_T)
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


    #ifdef NONTHERMAL_E_DIST
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
            if (i >=0  && i < N_PH_E && j >= 0 && j < N_GAMMA)
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
    #ifdef NONTHERMAL_E_DIST

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
