//
// Created by Tyler Parsotan on 11/18/25.
// This file holds the functions to create/read in the hot cross section lookup table
// This is calculated using the equations outlined in:
// Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009, Astrophysical Journal Supplement, 184, 387 &
// Tomohisa Kawashima et al 2023 ApJ 949 101 &
// Canfield E., Howard W. M. and Liang E. P. 1987 ApJ 323 565
//

#include "mcrat.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

// define the extent of the tabulated fluid photon energies normalized by electron rest mass
#define LOG_PH_E_MIN -12.0
#define LOG_PH_E_MAX -6.0
#define N_PH_E 220

// define the extent of the tabulated fluid temperatures normalized by electron rest mass
// for the thermal hot cross section calculation
#define LOG_T_MIN -4.0
#define LOG_T_MAX 4.0
#define N_T 80

//define the thermal modified cross section table filename
#define HOT_THERMAL_X_SECTION_FILE	"thermal_hot_x_section.dat"

//define the nonthermal modified cross section table filename
#define HOT_NONTHERMAL_X_SECTION_FILE	"nonthermal_hot_x_section.dat"

//helper struct to evaluate the double integral
struct double_integral_params { double norm_ph_comv; double theta };

double table[N_PH_E + 1][N_T + 1];

void initalizeHotCrossSection(gsl_rng *rand, FILE *fPtr)
{

    int i,j;
    double dt=(LOG_T_MAX-LOG_T_MIN)/N_T, dph_e=(LOG_PH_E_MAX-LOG_PH_E_MIN)/N_PH_E;
    double comv_ph_e, theta;
    FILE *fp;

    for (i = 0; i <= N_PH_E; i++)
    {
        for (j = 0; j <= N_T; j++)
        {
            comv_ph_e = pow(10., LOG_PH_E_MIN + i * dph_e);
            theta = pow(10., LOG_T_MIN + j * dt);
            table[i][j] = log10(calculateTotalThermalCrossSection(comv_ph_e, theta, rand, fPtr));
            if (isnan(table[i][j]))
            {
                fprintf(stdout, "%d %d %g %g\n", i, j, comv_ph_e, theta);
                exit(0);
            }
        }
    }

    fp = fopen(HOT_THERMAL_X_SECTION_FILE, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "couldn't write to file\n");
        exit(0);
    }
    for (i = 0; i <= N_PH_E; i++)
    {
        for (j = 0; j <= N_T; j++)
        {
            comv_ph_e = pow(10., LOG_PH_E_MIN + i * dph_e);
            theta = pow(10., LOG_T_MIN + j * dt);
            fprintf(fp, "%d %d %g %g %15.10g\n", i, j, comv_ph_e, theta, table[i][j]);
        }
    }
    fclose(fp);
    fprintf(stderr, "done.\n\n");

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
    double xl[2] = { -1, 1 };
    double xu[2] = { 1, 1. + 12 * theta };
    struct double_integral_params params = {ph_comv, theta };

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
    struct double_integral_params * fp = (struct double_integral_params *)p;
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

