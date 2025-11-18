//
// Created by Tyler Parsotan on 11/18/25.
// This file holds the functions to create/read in the hot cross section lookup table
// This is calculated using the equations outlined in:
// Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009, Astrophysical Journal Supplement, 184, 387 &
// Tomohisa Kawashima et al 2023 ApJ 949 101
//

double boostedCrossSection(double norm_ph_comv, double mu, double gamma)
{
    /*
        Calculates the KN cross section in the electron restframe
        This takes the photon's comobetaing energy (in the fluid frame) normalized by the electron rest mass, norm_ph_comv,
        the cosine of the angle between the photon and the electron in the electron rest frame, mu,
        and the lorentz factor of the electron in the fluid frame, gamma
    */
    double norm_ph_e=0, result=0, beta=0;

    /* energy in electron rest frame */
    beta = sqrt(gamma * gamma - 1.) / gamma;
    norm_ph_e = norm_ph_comv * gamma * (1. - mu * beta);

    result = hc_klein_nishina(norm_ph_e) * (1. - mu * beta);

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

