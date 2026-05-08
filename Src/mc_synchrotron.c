/*
 * mc_synchrotron.c
 * ================
 * Pitch-angle averaged synchrotron photon emission and SSA for non-thermal
 * electron distributions (single or broken power law).
 *
 * See mc_synchrotron.h for the physics overview and equation references.
 */

#include "mcrat.h"

/*
 * File-scope SynchUniversalTables instance.
 * Populated once by initSynchTables at startup (called from mcrat.c).
 * All functions in this file access it directly — no pointer passing needed.
 * Pattern mirrors global_interp_thermal_data in hot_x_section.c.
 */
static SynchUniversalTables synch_tables;
static int synch_tables_initialized = 0;

/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 1 — BESSEL INTEGRAL: F(x)                                          */
/* RAIKOU Eq. B13;  G&S91 Eq. 2;  R&L79 Eq. 6.31                             */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * bessel_K53_integrand
 * --------------------
 * Returns K_{5/3}(xi), the modified Bessel function of the second kind of
 * order 5/3, evaluated at xi. This is the integrand appearing in
 *
 *   F(x) = x * integral_x^inf K_{5/3}(xi) dxi         [RAIKOU Eq. B13]
 *
 * GSL's gsl_sf_bessel_Knu handles arbitrary real orders correctly.
 */
static double bessel_K53_integrand(double xi, void *params)
{
    (void)params;
    return gsl_sf_bessel_Knu(5.0/3.0, xi);
}


/*
 * computeFatX
 * -----------
 * Compute the synchrotron spectral function
 *
 *   F(x) = x * integral_x^inf K_{5/3}(xi) dxi         [RAIKOU Eq. B13]
 *
 * at a single value of the dimensionless frequency ratio
 *
 *   x = nu_f / (gamma_e^2 * nu_c)
 *
 * F(x) appears as the kernel of both the emission spectral integral G(x)
 * [RAIKOU Eq. B12] and the absorption spectral integral G_a(x)
 * [RAIKOU Eq. C3]. It is positive for all x > 0 and decays exponentially
 * for x >> 1 (R&L79 Eq. 6.33b).
 *
 * Parameters
 * ----------
 * x   : dimensionless frequency ratio, x > 0
 * ws  : pre-allocated GSL quadrature workspace (size >= 1000)
 * fPtr: log file
 *
 * Returns F(x) >= 0; returns 0.0 for x >= SYNCH_X_MAX.
 */
static double computeFatX(double x,
                           gsl_integration_workspace *ws,
                           FILE *fPtr)
{
    if (x >= SYNCH_X_MAX) return 0.0;

    gsl_function F;
    F.function = bessel_K53_integrand;
    F.params   = NULL;

    double result = 0.0, abserr = 0.0;
    int status = gsl_integration_qagiu(&F, x,
                                        1e-14, 1e-10,
                                        1000, ws,
                                        &result, &abserr);

    if (status != GSL_SUCCESS && x < 40.0)
    {
        fprintf(fPtr,
                ">> [computeFatX] WARNING: QAGIU status '%s' at "
                "x = %.4e (abserr = %.2e). "
                "Consider increasing SYNCH_N_X.\n",
                gsl_strerror(status), x, abserr);
        fflush(fPtr);
    }

    return x * result;
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 2 — ABSORPTION SPECTRAL INTEGRAL: G_a(x)                           */
/* RAIKOU Eq. C3                                                               */
/* ═══════════════════════════════════════════════════════════════════════════ */

/* Params struct for the G_a integrand */
struct Ga_integrand_params
{
    double            p;
    gsl_spline       *F_spline;
    gsl_interp_accel *F_acc;
};

/*
 * Ga_integrand_fn
 * ---------------
 * Evaluates the integrand z^{(p-2)/2} * F(z) for the absorption spectral
 * integral G_a(x; p) = integral_x^inf z^{(p-2)/2} F(z) dz [RAIKOU Eq. C3].
 */
static double Ga_integrand_fn(double z, void *vp)
{
    struct Ga_integrand_params *par = (struct Ga_integrand_params *)vp;
    double Fz = gsl_spline_eval(par->F_spline, z, par->F_acc);
    return pow(z, 0.5*(par->p - 2.0)) * Fz;
}


/*
 * computeGaAtX
 * ------------
 * Compute the absorption spectral integral
 *
 *   G_a(x; p) = integral_x^inf z^{(p-2)/2} F(z) dz    [RAIKOU Eq. C3]
 *
 * at a single value of x and power-law index p.
 *
 * This quantity depends only on x and p, not on B, nu_f, or n_e individually.
 * A 1D table of G_a(x) at fixed p therefore captures the entire spectral
 * shape of alpha_{nu_f}^(f), with B and n_e entering analytically through
 * the RAIKOU Eq. C2 prefactor evaluated at runtime.
 *
 * Parameters
 * ----------
 * x      : lower integration limit (dimensionless frequency ratio)
 * p      : power-law index of the electron distribution
 * F_spl  : pre-built spline of F(z) on [SYNCH_X_MIN, SYNCH_X_MAX]
 * F_acc  : GSL accelerator for F_spl
 * ws     : pre-allocated GSL quadrature workspace
 * fPtr   : log file
 *
 * Returns G_a(x; p) >= 0; returns 0.0 for x >= SYNCH_X_MAX.
 */
static double computeGaAtX(double x,
                             double p,
                             gsl_spline               *F_spl,
                             gsl_interp_accel         *F_acc,
                             gsl_integration_workspace *ws,
                             FILE *fPtr)
{
    if (x >= SYNCH_X_MAX) return 0.0;

    struct Ga_integrand_params par = { p, F_spl, F_acc };
    gsl_function G;
    G.function = Ga_integrand_fn;
    G.params   = &par;

    double result = 0.0, abserr = 0.0;
    double upper  = SYNCH_X_MAX;

    int status = gsl_integration_qags(&G, x, upper,
                                       1e-14, 1e-8,
                                       1000, ws,
                                       &result, &abserr);

    if (status != GSL_SUCCESS)
    {
        fprintf(fPtr,
                ">> [computeGaAtX] WARNING: QAGS status '%s' at "
                "x = %.4e, p = %.3f (abserr = %.2e).\n",
                gsl_strerror(status), x, p, abserr);
        fflush(fPtr);
    }

    return (result > 0.0) ? result : 0.0;
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 3 — UNIVERSAL TABLE INITIALISATION AND TEARDOWN                    */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * initSynchTables
 * ---------------
 * Build all x-dependent spectral tables. Must be called once before any
 * other function in this file.
 *
 * (1) F(x) grid  [RAIKOU Eq. B13]
 *     F_arr[i] = computeFatX(x_arr[i]) for i = 0..SYNCH_N_X-1
 *     on a log-spaced grid over [SYNCH_X_MIN, SYNCH_X_MAX].
 *
 * (2) G_a(x) grid  [RAIKOU Eq. C3]
 *     For POWERLAW:      Ga_arr[i]    = computeGaAtX(x_arr[i], POWERLAW_INDEX)
 *     For BROKENPOWERLAW: Ga_arr_p1[i] = computeGaAtX(x_arr[i], POWERLAW_INDEX_1)
 *                         Ga_arr_p2[i] = computeGaAtX(x_arr[i], POWERLAW_INDEX_2)
 *     Both tables are needed to evaluate RAIKOU Eq. C4.
 *
 * (3) Inverse CDF of x ~ F(x)*x d(log x)  [G&S91 Sec. 2]
 *     Allows direct sampling of x without rejection.
 *
 * (4) Inverse CDF of alpha ~ sin^2(alpha)  [G&S91 Sec. 2; R&L79 Sec. 6.2]
 *     CDF(alpha) = (alpha - sin(2*alpha)/2) / pi.
 */
void initSynchTables(FILE *fPtr)
{
    int i;
    fprintf(fPtr,
            ">> [initSynchTables] Building universal spectral tables "
            "(SYNCH_N_X = %d)...\n", SYNCH_N_X);
    fflush(fPtr);

    /* ── Allocate arrays ──────────────────────────────────────────────────── */
    synch_tables.x_arr       = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.F_arr       = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.Ga_arr      = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.Ga_arr_p1   = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.Ga_arr_p2   = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_x_cdf_u    = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_x_cdf_logx = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_alpha_cdf_u     = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_alpha_cdf_alpha = (double *)malloc(SYNCH_N_X * sizeof(double));

    /* ── Log-spaced x grid ───────────────────────────────────────────────── */
    double log_x_min = log10(SYNCH_X_MIN);
    double log_x_max = log10(SYNCH_X_MAX);
    for (i = 0; i < SYNCH_N_X; i++)
    {
        double t = (double)i / (SYNCH_N_X - 1);
        synch_tables.x_arr[i] = pow(10.0, log_x_min + t*(log_x_max - log_x_min));
    }

    /* ── (1) F(x)  [RAIKOU Eq. B13] ─────────────────────────────────────── */
    fprintf(fPtr,
            ">> [initSynchTables] Computing F(x) "
            "[RAIKOU Eq. B13]...\n");
    fflush(fPtr);

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);

    for (i = 0; i < SYNCH_N_X; i++)
        synch_tables.F_arr[i] = computeFatX(synch_tables.x_arr[i], ws, fPtr);

    synch_tables.F_acc    = gsl_interp_accel_alloc();
    synch_tables.F_spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    gsl_spline_init(synch_tables.F_spline,
                    synch_tables.x_arr, synch_tables.F_arr, SYNCH_N_X);

    /* ── (2) G_a(x)  [RAIKOU Eq. C3] ────────────────────────────────────── */
    fprintf(fPtr,
            ">> [initSynchTables] Computing G_a(x) "
            "[RAIKOU Eq. C3]...\n");
    fflush(fPtr);

    synch_tables.Ga_acc    = gsl_interp_accel_alloc();
    synch_tables.Ga_acc_p1 = gsl_interp_accel_alloc();
    synch_tables.Ga_acc_p2 = gsl_interp_accel_alloc();

    #if NONTHERMAL_E_DIST == POWERLAW

        for (i = 0; i < SYNCH_N_X; i++)
            synch_tables.Ga_arr[i] = computeGaAtX(synch_tables.x_arr[i],
                                              POWERLAW_INDEX,
                                                   synch_tables.F_spline,
                                                   synch_tables.F_acc,
                                              ws, fPtr);

        /* p1 / p2 arrays are unused for a single power law; zero-fill */
        for (i = 0; i < SYNCH_N_X; i++)
            synch_tables.Ga_arr_p1[i] = synch_tables.Ga_arr_p2[i] = 0.0;

        synch_tables.Ga_spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
        gsl_spline_init(synch_tables.Ga_spline,
                        synch_tables.x_arr, synch_tables.Ga_arr, SYNCH_N_X);

        /* Allocate p1/p2 splines as empty so freeSynchTables can always free */
        synch_tables.Ga_spline_p1 = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
        synch_tables.Ga_spline_p2 = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        /*
         * For the broken power law (RAIKOU Eq. C4) we need two separate G_a
         * tables — one evaluated with p1 (for the low-energy segment
         * gamma_min <= gamma <= gamma_br) and one with p2 (for the high-energy
         * segment gamma_br < gamma <= gamma_max). Ga_arr is set as an alias
         * for Ga_arr_p1 for consistency with the single power-law path.
         */
        for (i = 0; i < SYNCH_N_X; i++)
        {
            synch_tables.Ga_arr_p1[i] = computeGaAtX(synch_tables.x_arr[i],
                                                  POWERLAW_INDEX_1,
                                                      synch_tables.F_spline,
                                                      synch_tables.F_acc,
                                                  ws, fPtr);
            synch_tables.Ga_arr_p2[i] = computeGaAtX(synch_tables.x_arr[i],
                                                  POWERLAW_INDEX_2,
                                                      synch_tables.F_spline,
                                                      synch_tables.F_acc,
                                                  ws, fPtr);
            synch_tables.Ga_arr[i] = synch_tables.Ga_arr_p1[i];
        }

        synch_tables.Ga_spline    = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
        synch_tables.Ga_spline_p1 = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
        synch_tables.Ga_spline_p2 = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);

        gsl_spline_init(synch_tables.Ga_spline,
                        synch_tables.x_arr, synch_tables.Ga_arr, SYNCH_N_X);
        gsl_spline_init(synch_tables.Ga_spline_p1,
                        synch_tables.x_arr, synch_tables.Ga_arr_p1, SYNCH_N_X);
        gsl_spline_init(synch_tables.Ga_spline_p2,
                        synch_tables.x_arr, synch_tables.Ga_arr_p2, SYNCH_N_X);

    #endif

    /* ── (3) Inverse CDF of x ~ F(x)*x d(log x)  [G&S91 Sec. 2] ───────── */
    /*
     * The probability that a photon from a single electron has dimensionless
     * frequency ratio in [x, x+dx] is proportional to F(x) dx (R&L79
     * Eq. 6.36). Rewriting in d(log x) measure gives weight F(x)*x per
     * unit log x. We normalise the cumulative sum of this weight to [0,1]
     * and store the inverse map u -> log10(x) as a linear spline.
     */
    {
        double *Fx_weight = (double *)malloc(SYNCH_N_X * sizeof(double));
        double  dlog_x    = (log_x_max - log_x_min) / (SYNCH_N_X - 1);
        double  cum       = 0.0;

        for (i = 0; i < SYNCH_N_X; i++)
            Fx_weight[i] = synch_tables.F_arr[i] * synch_tables.x_arr[i];

        synch_tables.inv_x_cdf_u[0]    = 0.0;
        synch_tables.inv_x_cdf_logx[0] = log10(synch_tables.x_arr[0]);
        for (i = 1; i < SYNCH_N_X; i++)
        {
            cum += 0.5*(Fx_weight[i] + Fx_weight[i-1]) * dlog_x;
            synch_tables.inv_x_cdf_u[i]    = cum;
            synch_tables.inv_x_cdf_logx[i] = log10(synch_tables.x_arr[i]);
        }
        for (i = 0; i < SYNCH_N_X; i++)
            synch_tables.inv_x_cdf_u[i] /= cum;
        synch_tables.inv_x_cdf_u[SYNCH_N_X - 1] = 1.0;

        free(Fx_weight);
    }

    synch_tables.inv_x_cdf_acc    = gsl_interp_accel_alloc();
    synch_tables.inv_x_cdf_spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    gsl_spline_init(synch_tables.inv_x_cdf_spline,
                    synch_tables.inv_x_cdf_u,
                    synch_tables.inv_x_cdf_logx,
                    SYNCH_N_X);

    /* ── (4) Inverse CDF of alpha ~ sin^2(alpha)  [G&S91 Sec. 2] ───────── */
    /*
     * The isotropic pitch-angle distribution is
     *   f(alpha) = (2/pi) * sin^2(alpha),  alpha in [0, pi]
     *
     * Its analytic CDF is:
     *   C(alpha) = (alpha - sin(2*alpha)/2) / pi
     *
     * We evaluate this on a uniform alpha grid and store the inverse map
     * u -> alpha as a linear spline.
     */
    for (i = 0; i < SYNCH_N_X; i++)
    {
        double alpha = M_PI * (double)i / (SYNCH_N_X - 1);
        synch_tables.inv_alpha_cdf_alpha[i] = alpha;
        synch_tables.inv_alpha_cdf_u[i]     = (alpha - 0.5*sin(2.0*alpha)) / M_PI;
    }
    synch_tables.inv_alpha_cdf_u[SYNCH_N_X - 1] = 1.0;

    synch_tables.inv_alpha_cdf_acc    = gsl_interp_accel_alloc();
    synch_tables.inv_alpha_cdf_spline = gsl_spline_alloc(gsl_interp_linear,
                                                      SYNCH_N_X);
    gsl_spline_init(synch_tables.inv_alpha_cdf_spline,
                    synch_tables.inv_alpha_cdf_u,
                    synch_tables.inv_alpha_cdf_alpha,
                    SYNCH_N_X);

    gsl_integration_workspace_free(ws);
    
    synch_tables_initialized = 1;

    fprintf(fPtr,
            ">> [initSynchTables] All universal tables ready.\n\n");
    fflush(fPtr);
}


/*
 * freeSynchTables
 * ---------------
 * Release all heap memory allocated by initSynchTables.
 * Must be called once at end of simulation.
 */
void freeSynchTables(FILE *fPtr)
{
    if (!synch_tables_initialized)
    {
        fprintf(fPtr,
                ">> [freeSynchTables] WARNING: called before initSynchTables "
                "or already freed. No-op.\n");
        fflush(fPtr);
        return;
    }

    free(synch_tables.x_arr);
    free(synch_tables.F_arr);
    free(synch_tables.Ga_arr);
    free(synch_tables.Ga_arr_p1);
    free(synch_tables.Ga_arr_p2);
    free(synch_tables.inv_x_cdf_u);
    free(synch_tables.inv_x_cdf_logx);
    free(synch_tables.inv_alpha_cdf_u);
    free(synch_tables.inv_alpha_cdf_alpha);

    gsl_spline_free(synch_tables.F_spline);
    gsl_spline_free(synch_tables.Ga_spline);
    gsl_spline_free(synch_tables.Ga_spline_p1);
    gsl_spline_free(synch_tables.Ga_spline_p2);
    gsl_spline_free(synch_tables.inv_x_cdf_spline);
    gsl_spline_free(synch_tables.inv_alpha_cdf_spline);

    gsl_interp_accel_free(synch_tables.F_acc);
    gsl_interp_accel_free(synch_tables.Ga_acc);
    gsl_interp_accel_free(synch_tables.Ga_acc_p1);
    gsl_interp_accel_free(synch_tables.Ga_acc_p2);
    gsl_interp_accel_free(synch_tables.inv_x_cdf_acc);
    gsl_interp_accel_free(synch_tables.inv_alpha_cdf_acc);
    
    synch_tables_initialized = 0;
}

/*
 * getSynchTables
 * --------------
 * Return a const pointer to the file-scope synch_tables for use by
 * functions in this translation unit. Exits if called before
 * initSynchTables.
 */
static const SynchUniversalTables *getSynchTables(FILE *fPtr)
{
    if (!synch_tables_initialized)
    {
        fprintf(fPtr,
                ">> [getSynchTables] FATAL: synch_tables accessed before "
                "initSynchTables was called.\n");
        fflush(fPtr);
        exit(1);
    }
    return &synch_tables;
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 4 — EMISSION MONTE CARLO SAMPLERS                                  */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * synchSampleX
 * ------------
 * Sample x = nu_f / (gamma_e^2 * nu_c) from the single-electron synchrotron
 * spectrum F(x) via the precomputed inverse CDF  [G&S91 Sec. 2; R&L79 Eq. 6.36].
 *
 * The photon frequency is then assembled as:
 *   nu_f = x * gamma_e^2 * nu_c(B, alpha)
 * where nu_c = 3 e B sin(alpha) / (4 pi me c)          [RAIKOU Eq. B7]
 *
 * Parameters
 * ----------
 * tables : initialised SynchUniversalTables
 * u      : uniform random variate in (0, 1)
 *
 * Returns x > 0.
 */
double synchSampleX(double u)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    
    double u_lo = tables->inv_x_cdf_u[0];
    double u_hi = tables->inv_x_cdf_u[SYNCH_N_X - 1];
    if (u < u_lo + 1e-12) u = u_lo + 1e-12;
    if (u > u_hi - 1e-12) u = u_hi - 1e-12;

    double log10_x = gsl_spline_eval(tables->inv_x_cdf_spline,
                                      u, tables->inv_x_cdf_acc);
    return pow(10.0, log10_x);
}


/*
 * synchSampleAlpha
 * ----------------
 * Sample the electron pitch angle alpha from the isotropic distribution
 *   f(alpha) = (2/pi) * sin^2(alpha)                   [G&S91 Sec. 2]
 *
 * alpha enters the photon frequency as nu_c ∝ sin(alpha) * gamma^2 * B
 *                                                        [RAIKOU Eq. B7]
 *
 * Parameters
 * ----------
 * tables : initialised SynchUniversalTables
 * u      : uniform random variate in (0, 1)
 *
 * Returns alpha in (0, pi).
 */
double synchSampleAlpha(double u)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    
    double u_lo = tables->inv_alpha_cdf_u[0];
    double u_hi = tables->inv_alpha_cdf_u[SYNCH_N_X - 1];
    if (u < u_lo + 1e-12) u = u_lo + 1e-12;
    if (u > u_hi - 1e-12) u = u_hi - 1e-12;

    return gsl_spline_eval(tables->inv_alpha_cdf_spline,
                            u, tables->inv_alpha_cdf_acc);
}


/*
 * synchSampleGammaEmission
 * ------------------------
 * Sample the electron Lorentz factor gamma_e for synchrotron PHOTON EMISSION
 * using the emission-weighted distribution N(gamma) * gamma^2.
 *
 * The gamma^2 weighting arises because the total synchrotron power emitted
 * by a single electron scales as gamma^2 (R&L79 Eq. 6.38). A photon drawn
 * at random from the full emission distribution was therefore emitted by an
 * electron with Lorentz factor distributed as N(gamma)*gamma^2, not N(gamma).
 *
 * This function must NOT be replaced by samplePowerLaw or
 * sampleBrokenPowerLawSubgroup from electron.c, which sample from N(gamma)
 * directly (correct for scattering, not for emission).
 *
 * SINGLE POWER LAW  N(gamma) ∝ gamma^{-p}  [RAIKOU Eq. A2]:
 *   emission weight ∝ gamma^{2-p} = gamma^q,  q = 2-p
 *   inverse CDF: gamma(u) = (gmin^q + u*(gmax^q - gmin^q))^{1/q}
 *   special case q=0 (p=2): gamma(u) = gmin * (gmax/gmin)^u
 *
 * BROKEN POWER LAW  [RAIKOU Eq. A4-A5]:
 *   split at gamma_break; each segment sampled with its share of the
 *   total emission power (integral of N(gamma)*gamma^2 dgamma).
 *
 * Parameters
 * ----------
 * rand : GSL Mersenne Twister RNG
 *
 * Returns gamma_e in [GAMMA_MIN, GAMMA_MAX].
 */
double synchSampleGammaEmission(gsl_rng *rand)
{
    double u = gsl_rng_uniform_pos(rand);

    #if NONTHERMAL_E_DIST == POWERLAW

        double p    = POWERLAW_INDEX;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;
        double q    = 2.0 - p;   /* exponent of emission-weighted distribution */

        if (fabs(q) < 1e-6)
            return gmin * pow(gmax / gmin, u);   /* p == 2: log-uniform */

        double gmin_q = pow(gmin, q);
        double gmax_q = pow(gmax, q);
        return pow(gmin_q + u*(gmax_q - gmin_q), 1.0/q);

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        double p1   = POWERLAW_INDEX_1;
        double p2   = POWERLAW_INDEX_2;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;
        double gbr  = GAMMA_BREAK;
        double q1   = 2.0 - p1;
        double q2   = 2.0 - p2;

        /*
         * Continuity factor at gamma_break: enforces N(gamma) continuous at
         * gamma_br [RAIKOU Eq. A4-A5], consistent with electron.c convention.
         */
        double C_cont = pow(gbr, p2 - p1);

        /*
         * Total emission power in each segment:
         *   I1 = integral_{gmin}^{gbr}  gamma^{-p1} * gamma^2 dgamma
         *   I2 = C_cont * integral_{gbr}^{gmax} gamma^{-p2} * gamma^2 dgamma
         */
        double I1, I2;

        if (fabs(q1) < 1e-6)
            I1 = log(gbr / gmin);
        else
            I1 = (pow(gbr, q1) - pow(gmin, q1)) / q1;

        if (fabs(q2) < 1e-6)
            I2 = C_cont * log(gmax / gbr);
        else
            I2 = C_cont * (pow(gmax, q2) - pow(gbr, q2)) / q2;

        double f1 = I1 / (I1 + I2);   /* fraction of emission from low segment */

        if (u < f1)
        {
            /* Sample from low-energy segment */
            double u1 = u / f1;
            if (fabs(q1) < 1e-6)
                return gmin * pow(gbr / gmin, u1);
            double gmin_q1 = pow(gmin, q1);
            double gbr_q1  = pow(gbr,  q1);
            return pow(gmin_q1 + u1*(gbr_q1 - gmin_q1), 1.0/q1);
        }
        else
        {
            /* Sample from high-energy segment */
            double u2 = (u - f1) / (1.0 - f1);
            if (fabs(q2) < 1e-6)
                return gbr * pow(gmax / gbr, u2);
            double gbr_q2  = pow(gbr,  q2);
            double gmax_q2 = pow(gmax, q2);
            return pow(gbr_q2 + u2*(gmax_q2 - gbr_q2), 1.0/q2);
        }

    #endif
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 5 — ABSORPTION COEFFICIENT: RAIKOU APPENDIX C                      */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * evalGa
 * ------
 * Evaluate G_a(x; p) from a precomputed spline, clamping x to the table
 * range. Returns 0 for x >= SYNCH_X_MAX where F(z) = 0.
 *
 * This is the core table lookup called by synchAlphaNu to evaluate the
 * spectral integral differences [G_a(x_max) - G_a(x_min)] that appear in
 * RAIKOU Eqs. C2 and C4.
 */
static inline double evalGa(double x,
                              gsl_spline       *spline,
                              gsl_interp_accel *acc)
{
    if (x >= SYNCH_X_MAX) return 0.0;
    if (x <= SYNCH_X_MIN) return gsl_spline_eval(spline, SYNCH_X_MIN, acc);
    return gsl_spline_eval(spline, x, acc);
}


/*
 * synchAlphaNu
 * ------------
 * Evaluate the SSA absorption coefficient alpha_{nu_f}^(f) [cm^{-1}] in
 * the fluid rest frame at comoving frequency nu_f [Hz], for a cell with
 * magnetic field B [G] and nonthermal electron number density n_e_nth
 * [cm^{-3}].
 *
 * SINGLE POWER LAW  [RAIKOU Eq. C2]:
 *
 *   alpha_{nu_f}^(f) =
 *       (p-1)(p+2) n_{e,nth} e^2 nu_c
 *       / (4 sqrt(3) me c (gamma_min^{1-p} - gamma_max^{1-p}))
 *       * (nu_f / nu_c)^{-(p+4)/2}
 *       * [G_a(x_max; p) - G_a(x_min; p)]
 *
 * where (pitch-angle averaged, sin(theta_B) -> 1):
 *   nu_c  = 3 e B / (4 pi me c)                         [RAIKOU Eq. B7]
 *   x_min = nu_f / (gamma_max^2 * nu_c)
 *   x_max = nu_f / (gamma_min^2 * nu_c)
 *   G_a(x; p) = integral_x^inf z^{(p-2)/2} F(z) dz     [RAIKOU Eq. C3]
 *
 * BROKEN POWER LAW  [RAIKOU Eq. C4]:
 *
 *   alpha_{nu_f}^(f) =
 *       A n_{e,nth} e^2 / (4 sqrt(3) me c nu_c)
 *       * { (p1+2) * (nu_f/nu_c)^{-(p1+4)/2}
 *             * [G_a(x_max; p1) - G_a(x_br; p1)]
 *         + (p2+2) * C_cont * (nu_f/nu_c)^{-(p2+4)/2}
 *             * [G_a(x_br; p2) - G_a(x_min; p2)] }
 *
 * where A is the broken power-law normalisation (RAIKOU Eq. B17,
 * implemented in electron.c as brokenPowerLawNorm), C_cont =
 * gamma_br^{p2-p1} is the continuity factor [RAIKOU Eq. A4-A5], and:
 *   x_br  = nu_f / (gamma_br^2  * nu_c)
 *   x_min = nu_f / (gamma_max^2 * nu_c)   [smallest x: largest gamma]
 *   x_max = nu_f / (gamma_min^2 * nu_c)   [largest  x: smallest gamma]
 *
 * The B dependence enters both through nu_c ∝ B in the prefactor AND
 * through the x arguments of G_a. Both contributions are evaluated exactly
 * here — no B-scaling approximation is made.
 *
 * Parameters
 * ----------
 * nu_f    : comoving photon frequency [Hz]
 * B       : magnetic field strength in the cell [G]
 * n_e_nth : nonthermal electron number density [cm^{-3}]
 * tables  : initialised SynchUniversalTables
 * fPtr    : log file
 *
 * Returns alpha_{nu_f}^(f) >= 0  [cm^{-1}].
 */
double synchAlphaNu(double nu_f,
                    double B,
                    double n_e_nth,
                    FILE *fPtr)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    
    if (nu_f <= 0.0 || B <= 0.0 || n_e_nth <= 0.0) return 0.0;

    /*
     * Pitch-angle averaged critical frequency  [RAIKOU Eq. B7]:
     *   nu_c = 3 e B / (4 pi me c)
     * The sin(theta_B) factor is unity after averaging over the isotropic
     * pitch-angle distribution f(alpha) = (2/pi) sin^2(alpha).
     */
    double nu_c = (3.0 * CHARGE_EL * B) / (4.0 * M_PI * M_EL * C_LIGHT);
    if (nu_c <= 0.0) return 0.0;

    #if NONTHERMAL_E_DIST == POWERLAW

        double p    = POWERLAW_INDEX;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;

        /*
         * Normalisation denominator from RAIKOU Eq. C2:
         *   gamma_min^{1-p} - gamma_max^{1-p}
         * For p=1 this is zero; guard and return 0.
         */
        if (fabs(p - 1.0) < 1e-6)
        {
            fprintf(fPtr,
                    ">> [synchAlphaNu] WARNING: p = 1, degenerate normalisation. "
                    "Returning 0.\n");
            fflush(fPtr);
            return 0.0;
        }
        double denom = pow(gmin, 1.0 - p) - pow(gmax, 1.0 - p);
        if (fabs(denom) < 1e-300) return 0.0;

        /*
         * Dimensionless frequency arguments  [x = nu_f / (gamma^2 nu_c)]:
         *   x_min corresponds to gamma_max (largest gamma -> smallest x)
         *   x_max corresponds to gamma_min (smallest gamma -> largest x)
         */
        double x_min = nu_f / (gmax * gmax * nu_c);
        double x_max = nu_f / (gmin * gmin * nu_c);

        /*
         * Spectral integral difference  [RAIKOU Eq. C3]:
         *   delta_Ga = G_a(x_max; p) - G_a(x_min; p)
         * Note: G_a is monotonically decreasing in x, so delta_Ga >= 0 when
         * x_max >= x_min, which holds because gmin <= gmax.
         */
        double Ga_xmax  = evalGa(x_max, tables->Ga_spline, tables->Ga_acc);
        double Ga_xmin  = evalGa(x_min, tables->Ga_spline, tables->Ga_acc);
        double delta_Ga = Ga_xmax - Ga_xmin;
        if (delta_Ga <= 0.0) return 0.0;

        /*
         * Full expression  [RAIKOU Eq. C2]:
         *
         *   alpha = (p-1)(p+2) n_e e^2 nu_c
         *           / (4 sqrt(3) me c denom)
         *           * (nu_f/nu_c)^{-(p+4)/2}
         *           * delta_Ga
         */
        double prefac  = ((p - 1.0) * (p + 2.0) * n_e_nth
                          * CHARGE_EL * CHARGE_EL * nu_c)
                       / (4.0 * sqrt(3.0) * M_EL * C_LIGHT * denom);

        double spectral = pow(nu_f / nu_c, -0.5*(p + 4.0)) * delta_Ga;

        return prefac * spectral;

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        double p1   = POWERLAW_INDEX_1;
        double p2   = POWERLAW_INDEX_2;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;
        double gbr  = GAMMA_BREAK;

        /*
         * Normalisation constant A  [RAIKOU Eq. B17].
         * brokenPowerLawNorm is implemented in electron.c and returns the
         * constant A such that integral_{gmin}^{gmax} N(gamma) dgamma = 1
         * when N is expressed per unit n_e_nth.
         */
        double A = brokenPowerLawNorm(p1, p2, gmin, gmax, gbr);
        if (A <= 0.0) return 0.0;

        /*
         * Continuity factor gamma_br^{p2-p1} that enforces N(gamma) continuous
         * at gamma_br  [RAIKOU Eq. A4-A5].
         */
        double C_cont = pow(gbr, p2 - p1);

        /*
         * Dimensionless frequency arguments for both segments:
         *   x_max = nu_f / (gmin^2 * nu_c)   lower limit of low-segment integral
         *   x_br  = nu_f / (gbr^2  * nu_c)   boundary between segments
         *   x_min = nu_f / (gmax^2 * nu_c)   lower limit of high-segment integral
         */
        double x_max = nu_f / (gmin * gmin * nu_c);
        double x_br  = nu_f / (gbr  * gbr  * nu_c);
        double x_min = nu_f / (gmax * gmax * nu_c);

        /*
         * Spectral integral differences  [RAIKOU Eq. C3]:
         *
         * Low-energy segment (p1):
         *   delta_Ga_p1 = G_a(x_max; p1) - G_a(x_br; p1)
         *   spans gamma in [gmin, gbr], i.e. x in [x_br, x_max]
         *
         * High-energy segment (p2):
         *   delta_Ga_p2 = G_a(x_br; p2) - G_a(x_min; p2)
         *   spans gamma in [gbr, gmax], i.e. x in [x_min, x_br]
         */
        double Ga_xmax_p1 = evalGa(x_max, tables->Ga_spline_p1, tables->Ga_acc_p1);
        double Ga_xbr_p1  = evalGa(x_br,  tables->Ga_spline_p1, tables->Ga_acc_p1);
        double delta_Ga_p1 = Ga_xmax_p1 - Ga_xbr_p1;
        if (delta_Ga_p1 < 0.0) delta_Ga_p1 = 0.0;

        double Ga_xbr_p2  = evalGa(x_br,  tables->Ga_spline_p2, tables->Ga_acc_p2);
        double Ga_xmin_p2 = evalGa(x_min, tables->Ga_spline_p2, tables->Ga_acc_p2);
        double delta_Ga_p2 = Ga_xbr_p2 - Ga_xmin_p2;
        if (delta_Ga_p2 < 0.0) delta_Ga_p2 = 0.0;

        if (delta_Ga_p1 == 0.0 && delta_Ga_p2 == 0.0) return 0.0;

        /*
         * Full broken power-law expression  [RAIKOU Eq. C4]:
         *
         *   alpha = A n_{e,nth} e^2 / (4 sqrt(3) me c nu_c)
         *           * { (p1+2) * (nu_f/nu_c)^{-(p1+4)/2} * delta_Ga_p1
         *             + (p2+2) * C_cont * (nu_f/nu_c)^{-(p2+4)/2} * delta_Ga_p2 }
         */
        double common = (A * n_e_nth * CHARGE_EL * CHARGE_EL)
                      / (4.0 * sqrt(3.0) * M_EL * C_LIGHT * nu_c);

        double term1  = (p1 + 2.0)
                      * pow(nu_f / nu_c, -0.5*(p1 + 4.0))
                      * delta_Ga_p1;

        double term2  = (p2 + 2.0) * C_cont
                      * pow(nu_f / nu_c, -0.5*(p2 + 4.0))
                      * delta_Ga_p2;

        return common * (term1 + term2);

    #endif
}



/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 6 — SSA WEIGHT MODIFICATION                                         */
/* RAIKOU Eqs. 31, 37, 40                                                      */
/* ═══════════════════════════════════════════════════════════════════════════ */
/*
 * applySSAAbsorption
 * ------------------
 * Apply continuous SSA weight attenuation over a lab-frame step of length
 * dl [cm]:
 *
 *   w_new = w_old * exp(-ssa_abs_coeff * dl)     [RAIKOU Eq. 40]
 *
 * The ssa_abs_coeff was set by calculateOpticalDepth using the fluid-frame
 * opacity evaluated at the photon's comoving frequency and the local B and
 * n_e_nth [RAIKOU Eq. C2].
 *
 * Note: no Doppler factor (nu_f / nu_z) is applied here because dl is
 * already the lab-frame path length and the coeff was computed in the
 * fluid frame. The two frame corrections cancel for the weight update when
 * the medium is non-relativistic; for relativistic flows the full RAIKOU
 * Eq. 31 correction should be applied in calculateOpticalDepth instead.
 *
 * Parameters
 * ----------
 * ph : photon packet (modified in-place)
 * dl : lab-frame step length [cm]; must be >= 0
 */
void applyabsorption(struct photon *ph, double dl)
{
    #if SYNCHROTRON_SWITCH == ON
        if (ph == NULL)             return;
        if (dl <= 0.0)              return;
        if (ph->ssa_abs_coeff <= 0.0) return;

        double tau = ph->abs_optical_depth * dl;

        /* Guard against exp underflow for very large optical depths */
        if (tau > 700.0)
        {
            ph->weight = 0.0;
            return;
        }

        ph->weight *= exp(-tau);
    #else
        (void)ph;
        (void)dl;
    #endif
}

/*
 * calculateOpticalDepthSSA
 * ------------------------
 * Compute the SSA absorption coefficient alpha_{nu_f}^(f) [cm^{-1}] for
 * the photon's current cell and store it in ph->abs_optical_depth.
 *
 * Called from calculateOpticalDepth in optical_depth.c whenever
 * recalc_properties == 1, i.e. on the same trigger that refreshes the
 * scattering optical depth. Uses the file-scope synch_tables global
 * directly — no pointer passing needed, matching the hot_x_section.c
 * pattern.
 *
 * Only acts on SYNCH_PHOTON packets in cells where B > 0 and n_e_nth > 0.
 * For all other photon types abs_optical_depth is set to 0.0 so that
 * applyabsorption is a no-op.
 *
 * Parameters
 * ----------
 * ph         : photon packet (ph->abs_optical_depth modified in-place)
 * hydro_data : hydro frame (provides B and nonthermal_dens)
 * fPtr       : log file
 */
void calculateOpticalDepthSSA(struct photon          *ph,
                               struct hydro_dataframe *hydro_data,
                               FILE                   *fPtr)
{
    #if SYNCHROTRON_SWITCH == ON

    ph->abs_optical_depth = 0.0;

    if (ph->type != SYNCH_PHOTON)
        return;

    int ci = ph->nearest_block_index;
    if (ci < 0)
        return;

    double B_cell  = getMagneticFieldMagnitude(hydro_data, ci);
    double n_e_nth = (hydro_data->nonthermal_dens)[ci];

    if (B_cell <= 0.0 || n_e_nth <= 0.0)
        return;

    /*
     * Recover comoving frequency from comv_p0 [erg/c]:
     *   nu_f = comv_p0 * C_LIGHT / PL_CONST
     */
    double nu_f = ph->comv_p0 * C_LIGHT / PL_CONST;

    ph->abs_optical_depth = synchAlphaNu(nu_f, B_cell, n_e_nth, fPtr);

    #else
    (void)ph;
    (void)hydro_data;
    (void)fPtr;
    #endif
}
/*
 * applySynchSSAWeightModification
 * --------------------------------
 * Attenuate the weight of photon ph due to synchrotron self-absorption
 * along a comoving path of length dl [cm] through a cell.
 *
 * The optical depth increment for SSA absorption is  [RAIKOU Eq. 31]:
 *
 *   Delta tau_{nu_f}^(a) = (nu_f / nu_z) * alpha_{nu_f}^(f) * dl
 *
 * where nu_f is the photon frequency in the fluid rest frame, nu_z is the
 * photon frequency in the lab / ZAMO frame, and dl is the path length in
 * that same frame. The ratio (nu_f / nu_z) is the covariant frame correction
 * that transforms the fluid-frame opacity to the lab-frame optical depth
 * increment [RAIKOU Eq. 31].
 *
 * The weight is then modified by  [RAIKOU Eq. 40]:
 *
 *   w_new = w_old * exp(-Delta tau_{nu_f}^(a))
 *
 * This is the continuous-absorption treatment in which the photon packet
 * carries a reduced weight proportional to the survival probability, rather
 * than being stochastically destroyed. This is consistent with the approach
 * described in RAIKOU Sec. 5.2.
 *
 * Parameters
 * ----------
 * ph           : photon packet to modify in-place
 * dl           : lab-frame path length through the cell [cm]
 * B_cell       : magnetic field magnitude in the cell [G]
 * n_e_nth_cell : nonthermal electron number density in the cell [cm^{-3}]
 * tables       : initialised SynchUniversalTables
 * fPtr         : log file
 */
void applySynchSSAWeightModification(struct photon              *ph,
                                      double                      dl,
                                      double                      B_cell,
                                      double                      n_e_nth_cell,
                                      FILE                       *fPtr)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    
    if (dl <= 0.0 || B_cell <= 0.0 || n_e_nth_cell <= 0.0) return;

    /*
     * Comoving photon frequency nu_f [Hz]:
     *   ph->comv_p0 is the comoving 0-component of the photon 4-momentum
     *   in units of [erg/c], so nu_f = comv_p0 * c / h.
     */
    double nu_f = (ph->comv_p0 * C_LIGHT) / PL_CONST;

    /*
     * Lab-frame photon frequency nu_z [Hz]:
     *   ph->p0 is the lab-frame 0-momentum in units of [erg/c].
     */
    double nu_z = (ph->p0 * C_LIGHT) / PL_CONST;

    if (nu_f <= 0.0 || nu_z <= 0.0) return;

    /*
     * Absorption coefficient in the fluid rest frame [cm^{-1}]
     * following RAIKOU Eq. C2 (single power law) or C4 (broken power law).
     */
    double alpha_nu_f = synchAlphaNu(nu_f, B_cell, n_e_nth_cell,
                                    fPtr);

    /*
     * Optical depth increment  [RAIKOU Eq. 31]:
     *   Delta tau = (nu_f / nu_z) * alpha_{nu_f}^(f) * dl
     */
    double delta_tau = (nu_f / nu_z) * alpha_nu_f * dl;

    /*
     * Weight modification  [RAIKOU Eq. 40]:
     *   w_new = w_old * exp(-Delta tau)
     */
    ph->weight *= exp(-delta_tau);
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 7 — PER-CELL STRATIFIED FREQUENCY SAMPLER                          */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * synchNaturalNu  — unchanged from previously approved version.
 * See Section 7 documentation above.
 */
static double synchNaturalNu(double                      B_cell,
                               gsl_rng                    *rand)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    
    double x      = synchSampleX(tables, gsl_rng_uniform_pos(rand));
    double alpha  = synchSampleAlpha(tables, gsl_rng_uniform_pos(rand));
    double gamma_e;

    #if NONTHERMAL_E_DIST == POWERLAW
        gamma_e = samplePowerLaw(POWERLAW_INDEX, GAMMA_MIN, GAMMA_MAX,
                                  rand, NULL);
    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
        gamma_e = sampleBrokenPowerLawSubgroup(POWERLAW_INDEX_1, POWERLAW_INDEX_2,
                                                GAMMA_MIN, GAMMA_MAX, GAMMA_BREAK,
                                                rand, NULL);
    #endif

    double nu_c = (3.0 * CHARGE_EL * B_cell * sin(alpha) * gamma_e * gamma_e)
                / (4.0 * M_PI * M_EL * C_LIGHT);
    return x * nu_c;
}


/*
 * buildSynchCellStrata
 * --------------------
 * Build the stratified frequency sampling parameters for a single fluid cell
 * with magnetic field B_cell.
 *
 * Full algorithm description: see mc_synchrotron.h SynchCellStrata and the
 * previously approved Section 7 documentation. This version removes the
 * `continue` statement in the reference sample loop.
 *
 * The previous implementation used:
 *
 *   for (i = 0; i < SYNCH_N_REF; i++) {
 *       ...
 *       if (log_nu < log_nu_min || log_nu >= log_nu_max) continue;
 *       ...
 *   }
 *
 * This is replaced with an explicit `if (in_range)` block so that the loop
 * body that populates counts[] and hist[] only executes when the sampled
 * frequency actually falls within the cell's natural emission range. The
 * logic is identical; only the control flow structure changes.
 */
void buildSynchCellStrata(SynchCellStrata            *cs,
                           double                      B_cell,
                           gsl_rng                    *rand,
                           FILE                       *fPtr)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    
    int i, k;

    cs->B_cell = B_cell;

    /*
     * ── (1) Stratum frequency range for this cell ─────────────────────────
     *
     * Pitch-angle averaged critical frequency  [RAIKOU Eq. B7]:
     *   nu_c_cell = 3 e B_cell / (4 pi me c)
     *
     * The cell's natural emission spectrum spans
     *   [nu_c_cell * SYNCH_X_MIN * gamma_min^2,
     *    nu_c_cell * SYNCH_X_MAX * gamma_max^2]
     *
     * using the x-table range [SYNCH_X_MIN, SYNCH_X_MAX] as the full
     * extent of the F(x) kernel  [RAIKOU Eq. B13].
     */
    double nu_c_cell   = (3.0 * CHARGE_EL * B_cell)
                       / (4.0 * M_PI * M_EL * C_LIGHT);
    double nu_min_cell = nu_c_cell * SYNCH_X_MIN * GAMMA_MIN * GAMMA_MIN;
    double nu_max_cell = nu_c_cell * SYNCH_X_MAX * GAMMA_MAX * GAMMA_MAX;

    if (nu_min_cell <= 0.0 || nu_max_cell <= nu_min_cell)
    {
        /*
         * Degenerate frequency range — B_cell is effectively zero or the
         * gamma range produces a zero-width spectrum. Mark all strata as
         * empty and set the spline pointer to NULL so the caller can detect
         * this condition and apply the fallback frequency.
         */
        fprintf(fPtr,
                ">> [buildSynchCellStrata] WARNING: degenerate frequency "
                "range for B_cell = %.3e G "
                "(nu_min = %.3e, nu_max = %.3e). "
                "Setting all stratum probabilities to zero.\n",
                B_cell, nu_min_cell, nu_max_cell);
        fflush(fPtr);

        for (k = 0; k <= SYNCH_N_STRATA; k++)
            cs->strata_edges[k] = 0.0;
        for (k = 0; k < SYNCH_N_STRATA; k++)
            cs->stratum_probs[k] = 0.0;

        cs->inv_nu_cdf_spline = NULL;
        cs->inv_nu_cdf_acc    = NULL;
        return;
    }

    /* ── (2) Stratum boundaries ───────────────────────────────────────────── */
    double log_nu_min  = log10(nu_min_cell);
    double log_nu_max  = log10(nu_max_cell);
    double dlog_strata = (log_nu_max - log_nu_min) / SYNCH_N_STRATA;

    for (k = 0; k <= SYNCH_N_STRATA; k++)
        cs->strata_edges[k] = pow(10.0, log_nu_min + k * dlog_strata);

    /* ── (3) Reference sample ─────────────────────────────────────────────── */
    /*
     * Draw SYNCH_N_REF photon frequencies from the natural emission
     * distribution using this cell's B_cell via synchNaturalNu. For each
     * draw, check whether the frequency falls within [nu_min_cell,
     * nu_max_cell] and, only if it does, increment the coarse stratum
     * counter and the fine histogram bin.
     *
     * The in-range check previously used `continue` to skip out-of-range
     * samples. It is now expressed as an explicit `if (in_range)` block
     * that contains all histogram update logic. The semantic meaning is
     * identical: out-of-range samples contribute nothing to counts[] or
     * hist[] and only increment n_in_range when they are in range.
     */
    int    counts[SYNCH_N_STRATA];
    double hist[SYNCH_N_CDF_BINS];
    memset(counts, 0, sizeof(counts));
    memset(hist,   0, sizeof(hist));

    double dlog_fine  = (log_nu_max - log_nu_min) / SYNCH_N_CDF_BINS;
    int    n_in_range = 0;

    for (i = 0; i < SYNCH_N_REF; i++)
    {
        double nu     = synchNaturalNu(tables, B_cell, rand);
        double log_nu = log10(nu);

        /*
         * Only process this sample if it falls within the cell's frequency
         * range. This replaces the previous `continue` statement with an
         * explicit conditional block. All histogram update logic lives
         * inside this block so that it is clear exactly what happens to
         * in-range versus out-of-range samples.
         */
        if (log_nu >= log_nu_min && log_nu < log_nu_max)
        {
            n_in_range++;

            /* Coarse stratum bin index */
            int k_idx = (int)((log_nu - log_nu_min) / dlog_strata);
            if (k_idx < 0)               k_idx = 0;
            if (k_idx >= SYNCH_N_STRATA) k_idx = SYNCH_N_STRATA - 1;
            counts[k_idx]++;

            /* Fine histogram bin index for marginal CDF */
            int b = (int)((log_nu - log_nu_min) / dlog_fine);
            if (b < 0)                b = 0;
            if (b >= SYNCH_N_CDF_BINS) b = SYNCH_N_CDF_BINS - 1;
            hist[b] += 1.0;
        }
    }

    /* Stratum probabilities p_k = n_k / N_in_range  [see class docstring] */
    for (k = 0; k < SYNCH_N_STRATA; k++)
    {
        cs->stratum_probs[k] = (n_in_range > 0)
                               ? (double)counts[k] / n_in_range
                               : 0.0;
    }

    /* ── (4) Marginal CDF and inverse spline ──────────────────────────────── */
    cs->cdf_log_nu_edges[0] = log_nu_min;
    cs->cdf_log_nu_vals[0]  = 0.0;

    double cum = 0.0;
    for (i = 0; i < SYNCH_N_CDF_BINS; i++)
    {
        cum += hist[i];
        cs->cdf_log_nu_edges[i + 1] = log_nu_min + (i + 1) * dlog_fine;
        cs->cdf_log_nu_vals [i + 1] = cum;
    }

    if (cum > 0.0)
    {
        for (i = 0; i <= SYNCH_N_CDF_BINS; i++)
            cs->cdf_log_nu_vals[i] /= cum;
    }
    cs->cdf_log_nu_vals[SYNCH_N_CDF_BINS] = 1.0;

    cs->inv_nu_cdf_acc    = gsl_interp_accel_alloc();
    cs->inv_nu_cdf_spline = gsl_spline_alloc(gsl_interp_linear,
                                               SYNCH_N_CDF_BINS + 1);
    gsl_spline_init(cs->inv_nu_cdf_spline,
                    cs->cdf_log_nu_vals,
                    cs->cdf_log_nu_edges,
                    SYNCH_N_CDF_BINS + 1);
}


/*
 * freeSynchCellStrata  — unchanged from previously approved version.
 */
void freeSynchCellStrata(SynchCellStrata *cs)
{
    if (cs->inv_nu_cdf_spline != NULL)
    {
        gsl_spline_free(cs->inv_nu_cdf_spline);
        cs->inv_nu_cdf_spline = NULL;
    }
    if (cs->inv_nu_cdf_acc != NULL)
    {
        gsl_interp_accel_free(cs->inv_nu_cdf_acc);
        cs->inv_nu_cdf_acc = NULL;
    }
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 8 — SINGLE-PHOTON FILL HELPER                                       */
/* (unchanged from previously approved version — reproduced in full)           */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * synchFillPhoton
 * ---------------
 * Populate all fields of a struct photon for synchrotron emission from
 * fluid cell cell_idx at comoving frequency nu_f_comv [Hz] with weight
 * ph_weight. See previous approved version for full documentation.
 *
 * Steps:
 *   (1) Isotropic comoving 4-momentum with energy h*nu_f_comv
 *   (2) Lorentz boost to lab frame — exact pattern from mc_cyclosynch.c
 *       (hydroVectorToCartesian -> negate -> lorentzBoost 'p')
 *   (3) Birth position via hydroCoordinateToMcratCoordinate with random
 *       offset within cell bounding box
 *   (4) Stokes s0=1, s1=s2=s3=0 (unpolarised)  [G&S91 Sec. 2]
 *   (5) type=CS_POOL_PHOTON, num_scatt=0, recalc_properties=1
 */
static void synchFillPhoton(struct photon          *ph,
                              int                     cell_idx,
                              double                  nu_f_comv,
                              double                  ph_weight,
                              struct hydro_dataframe *hydro_data,
                              gsl_rng                *rand,
                              FILE                   *fPtr)
{
    double p_comv[4], p_lab[4];
    double boost[3];
    double position_phi;
    double cartesian_pos[3];

    /* ── (1) Isotropic comoving direction ─────────────────────────────────── */
    double com_v_phi   = samplePhotonPhi(rand, fPtr);
    double com_v_theta = gsl_rng_uniform(rand) * M_PI;
    double p_mag       = (PL_CONST * nu_f_comv) / C_LIGHT;

    p_comv[0] = p_mag;
    p_comv[1] = p_mag * sin(com_v_theta) * cos(com_v_phi);
    p_comv[2] = p_mag * sin(com_v_theta) * sin(com_v_phi);
    p_comv[3] = p_mag * cos(com_v_theta);

    ph->comv_p0 = p_comv[0];
    ph->comv_p1 = p_comv[1];
    ph->comv_p2 = p_comv[2];
    ph->comv_p3 = p_comv[3];

    /* ── (2) Lorentz boost to lab frame ──────────────────────────────────── */
    /*
     * Exact pattern from mc_cyclosynch.c: populate boost[] with the
     * Cartesian fluid velocity, negate all three spatial components, then
     * call lorentzBoost with flag 'p' (4-momentum).
     */
    #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
        position_phi = gsl_rng_uniform(rand) * 2.0 * M_PI;
        hydroVectorToCartesian(boost,
                                (hydro_data->v0)[cell_idx],
                                (hydro_data->v1)[cell_idx],
                                (hydro_data->v2)[cell_idx],
                                (hydro_data->r0)[cell_idx],
                                (hydro_data->r1)[cell_idx],
                                position_phi);
    #else
        position_phi = 0.0;
        hydroVectorToCartesian(boost,
                                (hydro_data->v0)[cell_idx],
                                (hydro_data->v1)[cell_idx],
                                (hydro_data->v2)[cell_idx],
                                (hydro_data->r0)[cell_idx],
                                (hydro_data->r1)[cell_idx],
                                (hydro_data->r2)[cell_idx]);
    #endif
    boost[0] *= -1.0;
    boost[1] *= -1.0;
    boost[2] *= -1.0;

    lorentzBoost(boost, p_comv, p_lab, 'p', fPtr);

    ph->p0 = p_lab[0];
    ph->p1 = p_lab[1];
    ph->p2 = p_lab[2];
    ph->p3 = p_lab[3];

    /* ── (3) Birth position ──────────────────────────────────────────────── */
    double pos_rand0 = gsl_rng_uniform_pos(rand)
                     * (hydro_data->r0_size)[cell_idx]
                     - 0.5*(hydro_data->r0_size)[cell_idx];
    double pos_rand1 = gsl_rng_uniform_pos(rand)
                     * (hydro_data->r1_size)[cell_idx]
                     - 0.5*(hydro_data->r1_size)[cell_idx];

    #if DIMENSIONS == THREE
        double pos_rand2 = gsl_rng_uniform_pos(rand)
                         * (hydro_data->r2_size)[cell_idx]
                         - 0.5*(hydro_data->r2_size)[cell_idx];
        hydroCoordinateToMcratCoordinate(cartesian_pos,
                                          (hydro_data->r0)[cell_idx] + pos_rand0,
                                          (hydro_data->r1)[cell_idx] + pos_rand1,
                                          (hydro_data->r2)[cell_idx] + pos_rand2);
    #else
        hydroCoordinateToMcratCoordinate(cartesian_pos,
                                          (hydro_data->r0)[cell_idx] + pos_rand0,
                                          (hydro_data->r1)[cell_idx] + pos_rand1,
                                          position_phi);
    #endif

    ph->r0 = cartesian_pos[0];
    ph->r1 = cartesian_pos[1];
    ph->r2 = cartesian_pos[2];

    /* ── (4) Stokes: unpolarised  [G&S91 Sec. 2] ────────────────────────── */
    ph->s0 = 1.0;
    ph->s1 = 0.0;
    ph->s2 = 0.0;
    ph->s3 = 0.0;

    /* ── (5) Bookkeeping ─────────────────────────────────────────────────── */
    ph->type                = SYNCH_PHOTON;
    ph->num_scatt           = 0.0;
    ph->weight              = ph_weight;
    ph->nearest_block_index = cell_idx;
    ph->recalc_properties   = 1;
    ph->time_to_scatter     = 0.0;
    ph->total_optical_depth = 0.0;

    #if SCATTERING_BIAS_SWITCH != OFF
        {
            int idx;
            for (idx = 0; idx < 1 + N_GAMMA; idx++)
            {
                ph->optical_depths[idx]  = 0.0;
                ph->scattering_bias[idx] = 0.0;
            }
        }
    #endif
    
    #if SYNCHROTRON_SWITCH == ON
        ph->abs_optical_depth=0.0; //this holds the opacity for synchrotron absorption (multiply by path length to get optical depth when modifying the weight correspondent with synchrotron self absorption
    #endif

}




/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 9 — MAIN EMISSION FUNCTION: HYBRID PER-CELL PHYSICAL COUNT         */
/*             WITH PER-CELL STRATIFIED FREQUENCY SAMPLING                    */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * emitCellPackets
 * ---------------
 * Helper: emit all ph_count photon packets for a single active cell ci,
 * using the pre-built per-cell strata cs and placing results into
 * ph_emit[*idx_ptr ... *idx_ptr + ph_count - 1].
 *
 * This helper exists to keep photonEmitSynch readable: the per-cell
 * emission logic (strata validity check, stratum allocation, per-stratum
 * frequency sampling, SSA attenuation) is isolated here so that the main
 * function body is a clean loop over active cells.
 *
 * FALLBACK PATH
 * -------------
 * If cs->inv_nu_cdf_spline == NULL (degenerate B_cell) OR if no stratum
 * has p_k >= SYNCH_P_MIN, all ph_count packets are emitted at the fallback
 * frequency:
 *   nu_fallback = nu_c(B_cell) * GAMMA_MAX^2
 * with weight ph_weight_adjusted. This preserves ph_tot and photon list
 * integrity while logging the condition.
 *
 * NORMAL PATH
 * -----------
 * ph_count packets are divided across active strata using integer division
 * with remainder distributed to the lowest strata first, so that
 *   sum_k stratum_alloc[k] == ph_count  exactly.
 *
 * For each packet in stratum k:
 *   - nu_f drawn via inverse CDF within stratum  [RAIKOU Eq. B11]
 *   - weight = ph_weight_adjusted * p_k * SYNCH_N_STRATA
 *   - synchFillPhoton sets kinematics
 *   - applySynchSSAWeightModification attenuates at birth  [RAIKOU Eq. 40]
 *
 * Parameters
 * ----------
 * ci                 : hydro_data cell index
 * ph_count           : number of packets to emit from this cell
 * cs                 : per-cell strata (may have NULL spline)
 * ph_weight_adjusted : converged photon packet weight
 * ph_emit            : output photon array (must have space for ph_count)
 * idx_ptr            : pointer to the running write index into ph_emit
 * hydro_data         : fluid grid
 * tables             : SynchUniversalTables (for SSA)
 * rand               : GSL RNG
 * fPtr               : log file
 */
static void emitCellPackets(int                         ci,
                              int                         ph_count,
                              const SynchCellStrata      *cs,
                              double                      ph_weight_adjusted,
                              struct photon              *ph_emit,
                              int                        *idx_ptr,
                              struct hydro_dataframe     *hydro_data,
                              gsl_rng                    *rand,
                              FILE                       *fPtr)
{
    int    i, k;
    double B_cell  = getMagneticFieldMagnitude(hydro_data, ci);

    /*
     * Determine whether the normal stratified path is available.
     * Two conditions must both hold:
     *   (a) The strata were built successfully (spline is non-NULL).
     *   (b) At least one stratum has p_k >= SYNCH_P_MIN.
     * If either fails we use the fallback path.
     */
    int n_active_strata = 0;

    if (cs->inv_nu_cdf_spline != NULL)
    {
        for (k = 0; k < SYNCH_N_STRATA; k++)
        {
            if (cs->stratum_probs[k] >= SYNCH_P_MIN)
                n_active_strata++;
        }
    }

    if (n_active_strata == 0)
    {
        /*
         * FALLBACK PATH
         * -------------
         * Emit all ph_count packets at a single representative frequency:
         *   nu_fallback = nu_c(B_cell) * GAMMA_MAX^2
         * which is the critical frequency for the highest-energy electrons.
         * Weight is ph_weight_adjusted (no importance correction needed
         * since all packets are identical).
         *
         * This path is taken when:
         *   - B_cell is effectively zero (cs->inv_nu_cdf_spline == NULL), OR
         *   - the reference sample found no photons in any stratum (all
         *     natural emission lies outside [nu_min_cell, nu_max_cell]).
         */
        double nu_c_fallback = (3.0 * CHARGE_EL * B_cell)
                             / (4.0 * M_PI * M_EL * C_LIGHT);
        double nu_fallback   = nu_c_fallback * GAMMA_MAX * GAMMA_MAX;

        fprintf(fPtr,
                ">> [emitCellPackets] Cell %d: using fallback frequency "
                "nu = %.3e Hz (B_cell = %.3e G, n_active_strata = 0). "
                "Emitting %d packets.\n",
                ci, nu_fallback, B_cell, ph_count);
        fflush(fPtr);

        for (i = 0; i < ph_count; i++)
        {
            synchFillPhoton(&ph_emit[*idx_ptr], ci, nu_fallback,
                             ph_weight_adjusted, hydro_data, rand, fPtr);
            //applySynchSSAWeightModification(&ph_emit[*idx_ptr],
            //                                dl_birth, B_cell, n_e_nth,
            //                               tables, fPtr);
            (*idx_ptr)++;
        }
    }
    else
    {
        /*
         * NORMAL STRATIFIED PATH
         * ----------------------
         * Divide ph_count packets across active strata using integer
         * division; distribute any remainder (ph_count % n_active_strata)
         * one extra packet to the lowest-indexed active strata in order.
         * This guarantees sum_k stratum_alloc[k] == ph_count exactly, with
         * no packet lost to truncation.
         */
        int base_per_stratum = ph_count / n_active_strata;
        int remainder        = ph_count % n_active_strata;
        int remainder_given  = 0;

        int stratum_alloc[SYNCH_N_STRATA];
        for (k = 0; k < SYNCH_N_STRATA; k++)
        {
            if (cs->stratum_probs[k] >= SYNCH_P_MIN)
            {
                stratum_alloc[k] = base_per_stratum;
                if (remainder_given < remainder)
                {
                    stratum_alloc[k]++;
                    remainder_given++;
                }
            }
            else
            {
                stratum_alloc[k] = 0;
            }
        }

        for (k = 0; k < SYNCH_N_STRATA; k++)
        {
            /*
             * Only emit into stratum k if it received a non-zero allocation.
             * This replaces the previous `continue` with an explicit
             * conditional block that wraps all per-stratum emission logic.
             */
            if (stratum_alloc[k] > 0)
            {
                /*
                 * Importance weight for stratum k  [RAIKOU Eq. B11]:
                 *   w_k = ph_weight_adjusted * p_k * SYNCH_N_STRATA
                 *
                 * This corrects for the forced equal allocation across strata
                 * (each stratum gets 1/N_STRATA of ph_count packets, but the
                 * natural fraction is p_k). The weighted sum over all packets
                 * therefore recovers the correct RAIKOU Eq. B11 spectral shape.
                 */
                double w_stratum = ph_weight_adjusted
                                 * cs->stratum_probs[k]
                                 * (double)SYNCH_N_STRATA;

                /*
                 * CDF limits for conditional frequency sampling within
                 * stratum k. Scan the marginal CDF array to find the u
                 * values that bracket [log_nu_lo, log_nu_hi].
                 *
                 * The scan replaces the use of a separate forward-CDF spline
                 * and is O(SYNCH_N_CDF_BINS) per stratum per cell. Since
                 * this is called at most SYNCH_N_STRATA times per cell and
                 * SYNCH_N_CDF_BINS = 1000, the total cost per cell is at
                 * most 10000 comparisons — negligible compared to the
                 * SYNCH_N_REF reference draws in buildSynchCellStrata.
                 */
                double log_nu_lo  = log10(cs->strata_edges[k]);
                double log_nu_hi  = log10(cs->strata_edges[k + 1]);
                double u_lo       = cs->cdf_log_nu_vals[0];
                double u_hi       = cs->cdf_log_nu_vals[SYNCH_N_CDF_BINS];
                int    found_lo   = 0;
                int    found_hi   = 0;

                for (i = 0; i <= SYNCH_N_CDF_BINS; i++)
                {
                    if (found_lo == 0 &&
                        cs->cdf_log_nu_edges[i] >= log_nu_lo)
                    {
                        u_lo    = cs->cdf_log_nu_vals[i];
                        found_lo = 1;
                    }
                    if (found_hi == 0 &&
                        cs->cdf_log_nu_edges[i] >= log_nu_hi)
                    {
                        u_hi    = cs->cdf_log_nu_vals[i];
                        found_hi = 1;
                    }
                }

                /* Clamp to valid spline domain */
                if (u_lo < cs->cdf_log_nu_vals[0])
                    u_lo = cs->cdf_log_nu_vals[0];
                if (u_hi > cs->cdf_log_nu_vals[SYNCH_N_CDF_BINS])
                    u_hi = cs->cdf_log_nu_vals[SYNCH_N_CDF_BINS];

                /*
                 * Only emit packets for this stratum if the CDF interval
                 * has non-zero width. A zero-width interval means the
                 * stratum's frequency range has no support in the CDF,
                 * which can occur at the extreme tails of the spectrum
                 * where the reference sample has no counts. In this case
                 * we absorb these packets into the fallback frequency rather
                 * than losing them, preserving ph_count exactly.
                 *
                 * This replaces the previous `continue` on the degenerate
                 * CDF interval with an if/else that either samples normally
                 * or uses the fallback, so every allocated packet is emitted.
                 */
                if (u_hi > u_lo)
                {
                    for (i = 0; i < stratum_alloc[k]; i++)
                    {
                        /*
                         * Draw u ~ Uniform[u_lo, u_hi] and map through the
                         * inverse marginal CDF spline to get log10(nu_f)
                         * within stratum k. This is equivalent to drawing
                         * from RAIKOU Eq. B11 restricted to the stratum
                         * frequency interval.
                         */
                        double u_nu = u_lo
                                    + gsl_rng_uniform_pos(rand)
                                    * (u_hi - u_lo);

                        double nu_f = pow(10.0,
                                          gsl_spline_eval(
                                              cs->inv_nu_cdf_spline,
                                              u_nu,
                                              cs->inv_nu_cdf_acc));

                        /*
                         * Hard clamp to stratum bounds to guard against
                         * spline overshoot at the piecewise-linear CDF
                         * boundaries.
                         */
                        if (nu_f < cs->strata_edges[k])
                            nu_f = cs->strata_edges[k];
                        if (nu_f > cs->strata_edges[k + 1])
                            nu_f = cs->strata_edges[k + 1];

                        /*
                         * Fill all photon fields: comoving 4-momentum,
                         * Lorentz boost to lab frame, birth position, Stokes
                         * parameters (s0=1, s1=s2=s3=0), and bookkeeping
                         * (type=CS_POOL_PHOTON, recalc_properties=1).
                         */
                        synchFillPhoton(&ph_emit[*idx_ptr], ci, nu_f,
                                         w_stratum, hydro_data, rand, fPtr);

                        /*
                         * SSA pre-attenuation at birth  [RAIKOU Eqs. 31, 40]:
                         *   w_new = w_old * exp(-(nu_f/nu_z) * alpha * dl)
                         * dl = 0.5 * r0_size[ci] is the mean path length
                         * from a uniformly distributed birth point to the
                         * nearest cell boundary.
                         */
                        //applySynchSSAWeightModification(&ph_emit[*idx_ptr],
                        //                                dl_birth,
                        //                                B_cell, n_e_nth,
                        //                                tables, fPtr);
                        (*idx_ptr)++;
                    }
                }
                else
                {
                    /*
                     * Degenerate CDF interval: emit stratum_alloc[k] packets
                     * at the fallback frequency with ph_weight_adjusted so
                     * that ph_count is preserved exactly. Log the condition
                     * so it can be diagnosed if the spectrum looks wrong.
                     */
                    double nu_c_fallback = (3.0 * CHARGE_EL * B_cell)
                                         / (4.0 * M_PI * M_EL * C_LIGHT);
                    double nu_fallback   = nu_c_fallback
                                        * GAMMA_MAX * GAMMA_MAX;

                    fprintf(fPtr,
                            ">> [emitCellPackets] Cell %d stratum %d: "
                            "degenerate CDF interval [u_lo=%.4f, u_hi=%.4f]. "
                            "Emitting %d packets at fallback nu = %.3e Hz.\n",
                            ci, k, u_lo, u_hi, stratum_alloc[k], nu_fallback);
                    fflush(fPtr);

                    for (i = 0; i < stratum_alloc[k]; i++)
                    {
                        synchFillPhoton(&ph_emit[*idx_ptr], ci, nu_fallback,
                                         ph_weight_adjusted,
                                         hydro_data, rand, fPtr);
                        
                        (*idx_ptr)++;
                    }
                }
            }
        }
    }
}


/*
 * photonEmitSynch
 * ---------------
 * Emit non-thermal synchrotron photon packets from all active fluid cells
 * within a thin shell at the injection radius.
 *
 * Full physics documentation: see previously approved Section 9 docstring.
 * This version:
 *   (a) Removes all `continue` statements.
 *   (b) Delegates all per-cell packet emission to emitCellPackets, which
 *       itself contains no `continue` statements.
 *   (c) Uses explicit `if (ph_dens[j] > 0)` guards wherever the previous
 *       version used `continue` to skip zero-packet cells.
 *
 * STRUCTURE OVERVIEW
 * ------------------
 * Step 1 : Compute shell rmin / rmax via calcCyclosynchRLimits.
 * Step 2 : First pass over hydro grid — count active cells (block_cnt).
 * Step 3 : Allocate per-cell arrays ph_dens[], cell_index[], W_cell[].
 * Step 4 : Second pass — compute emission weights W_i = n_e * B^2 * V.
 * Step 5 : Weight-tuning loop — adjust ph_weight_adjusted until
 *           min_photons <= sum lambda_i <= max_photons.
 * Step 6 : Poisson draw — n_cell_i ~ Poisson(lambda_i) for each cell.
 * Step 7 : Allocate ph_emit[ph_tot] and initialise all entries.
 * Step 8 : Per-cell emission loop — for each cell j with ph_dens[j] > 0,
 *           build per-cell strata, call emitCellPackets.
 * Step 9 : Add all emitted photons to the photon list in a single batch call.
 * Step 10: Free all heap allocations.
 */
int photonEmitSynch(struct photonList          *photon_list,
                    double                      r_inj,
                    double                      ph_weight,
                    int                         min_photons,
                    int                         max_photons,
                    double                      theta_min,
                    double                      theta_max,
                    struct hydro_dataframe     *hydro_data,
                    gsl_rng                    *rand,
                    FILE                       *fPtr)
{
    int    i, j;
    int    n_emitted = 0;
    int    block_cnt = 0;

    double r_grid_innercorner    = 0.0;
    double r_grid_outercorner    = 0.0;
    double theta_grid_innercorner = 0.0;
    double theta_grid_outercorner = 0.0;

    /* ── Step 1: Shell boundaries  [mirrors photonEmitCyclosynch] ────────────
     *
     * rmin = r_inj - c/fps,  rmax = r_inj + c/fps.
     * The shell thickness equals the distance light travels between
     * consecutive hydro frames, capturing exactly the region of newly
     * emitting gas at each injection epoch.
     */
    double rmin = calcCyclosynchRLimits(hydro_data->scatt_frame_number,
                                         hydro_data->inj_frame_number,
                                         hydro_data->fps, r_inj, "min");
    double rmax = calcCyclosynchRLimits(hydro_data->scatt_frame_number,
                                             hydro_data->inj_frame_number,
                                             hydro_data->fps, r_inj, "max");

    fprintf(fPtr,
            ">> [photonEmitSynch] Shell: rmin = %.3e cm, "
            "rmax = %.3e cm, theta = [%.3f, %.3f] rad\n",
            rmin, rmax, theta_min, theta_max);
    fflush(fPtr);

    /* ── Step 2: First pass — count active cells (block_cnt) ─────────────────
     *
     * A cell is active if ALL of the following hold:
     *   (a) nonthermal_dens[i] > 0  (has a non-thermal electron population)
     *   (b) Its bounding-box corners overlap the shell [rmin, rmax] in r
     *   (c) Its bounding-box corners overlap [theta_min, theta_max] in theta
     *
     * The corner test uses hydroCoordinateToSpherical on the min and max
     * corners of the cell bounding box, identical to photonEmitCyclosynch.
     * Cells that only partially overlap the shell boundary are included,
     * which avoids a systematic undercount at the shell edges.
     *
     * block_cnt is used to allocate ph_dens[], cell_index[], and W_cell[]
     * to exactly the right size before the second pass.
     *
     * The previous implementation used `continue` inside the loop to skip
     * cells that fail the nonthermal_dens or geometry test. Here the
     * active-cell test is expressed as a single compound `if` whose body
     * increments block_cnt, so the loop has one clear execution path per
     * cell with no early-exit jumps.
     */
    for (i = 0; i < hydro_data->num_elements; i++)
    {
        int has_nonthermal = ((hydro_data->nonthermal_dens)[i] > 0.0);

        if (has_nonthermal)
        {
            #if DIMENSIONS == THREE
                hydroCoordinateToSpherical(&r_grid_innercorner,
                                            &theta_grid_innercorner,
                                            fabs((hydro_data->r0)[i])
                                                - 0.5*(hydro_data->r0_size)[i],
                                            fabs((hydro_data->r1)[i])
                                                - 0.5*(hydro_data->r1_size)[i],
                                            fabs((hydro_data->r2)[i])
                                                - 0.5*(hydro_data->r2_size)[i]);
                hydroCoordinateToSpherical(&r_grid_outercorner,
                                            &theta_grid_outercorner,
                                            fabs((hydro_data->r0)[i])
                                                + 0.5*(hydro_data->r0_size)[i],
                                            fabs((hydro_data->r1)[i])
                                                + 0.5*(hydro_data->r1_size)[i],
                                            fabs((hydro_data->r2)[i])
                                                + 0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner,
                                            &theta_grid_innercorner,
                                            (hydro_data->r0)[i]
                                                - 0.5*(hydro_data->r0_size)[i],
                                            (hydro_data->r1)[i]
                                                - 0.5*(hydro_data->r1_size)[i],
                                            0.0);
                hydroCoordinateToSpherical(&r_grid_outercorner,
                                            &theta_grid_outercorner,
                                            (hydro_data->r0)[i]
                                                + 0.5*(hydro_data->r0_size)[i],
                                            (hydro_data->r1)[i]
                                                + 0.5*(hydro_data->r1_size)[i],
                                            0.0);
            #endif
            int in_shell = (rmin <= r_grid_outercorner)
                        && (r_grid_innercorner  <  rmax)
                        && (theta_grid_outercorner >= theta_min)
                        && (theta_grid_innercorner  <  theta_max);

            if (in_shell)
                block_cnt++;
        }
    }

    fprintf(fPtr,
            ">> [photonEmitSynch] Found %d active cells in the shell.\n",
            block_cnt);
    fflush(fPtr);

    if (block_cnt == 0)
    {
        fprintf(fPtr,
                ">> [photonEmitSynch] WARNING: no active cells with "
                "nonthermal electrons in the injection shell. "
                "No photons emitted.\n\n");
        fflush(fPtr);
        return 0;
    }

    /* ── Step 3: Allocate per-cell arrays ────────────────────────────────────
     *
     * ph_dens[j]    : Poisson-drawn packet count for active cell j
     * cell_index[j] : hydro_data array index of active cell j
     * W_cell[j]     : unnormalised emission weight W_i = n_e_nth * B^2 * V
     *
     * All three arrays are indexed by j in [0, block_cnt), where j
     * increments only for cells that pass the active-cell test in Step 4.
     */
    int    *ph_dens    = (int    *)malloc(block_cnt * sizeof(int));
    int    *cell_index = (int    *)malloc(block_cnt * sizeof(int));
    double *W_cell     = (double *)malloc(block_cnt * sizeof(double));

    /* ── Step 4: Second pass — record cell indices and emission weights ───────
     *
     * W_i = n_{e,nth,i} * B_i^2 * V_i
     *
     * This is the unnormalised total synchrotron emission power of cell i,
     * obtained by integrating RAIKOU Eq. B11 over nu_f and gamma_e.
     * The B^2 factor arises from:
     *   integral_0^inf j_{nu_f}^(f) d nu_f
     *     ∝ n_{e,nth} * B^2 * integral N_hat(gamma) gamma^2 dgamma
     *                                      [R&L79 Eq. 6.38; G&S91 Eq. 4]
     * The gamma integral is cell-independent (depends only on GAMMA_MIN,
     * GAMMA_MAX, and POWERLAW_INDEX, all compile-time constants) and
     * factors out of the relative weight comparison between cells.
     *
     * hydroElementVolume (geometry.c) returns the proper coordinate-system
     * volume V_i accounting for the grid geometry and dimensionality.
     *
     * The same compound `if (has_nonthermal && in_shell)` structure as
     * Step 2 is used, with `j` incrementing only for active cells to keep
     * cell_index[], W_cell[], and ph_dens[] aligned.
     */
    j = 0;
    for (i = 0; i < hydro_data->num_elements; i++)
    {
        int has_nonthermal = ((hydro_data->nonthermal_dens)[i] > 0.0);

        if (has_nonthermal)
        {
            #if DIMENSIONS == THREE
                hydroCoordinateToSpherical(&r_grid_innercorner,
                                            &theta_grid_innercorner,
                                            fabs((hydro_data->r0)[i])
                                                - 0.5*(hydro_data->r0_size)[i],
                                            fabs((hydro_data->r1)[i])
                                                - 0.5*(hydro_data->r1_size)[i],
                                            fabs((hydro_data->r2)[i])
                                                - 0.5*(hydro_data->r2_size)[i]);
                hydroCoordinateToSpherical(&r_grid_outercorner,
                                            &theta_grid_outercorner,
                                            fabs((hydro_data->r0)[i])
                                                + 0.5*(hydro_data->r0_size)[i],
                                            fabs((hydro_data->r1)[i])
                                                + 0.5*(hydro_data->r1_size)[i],
                                            fabs((hydro_data->r2)[i])
                                                + 0.5*(hydro_data->r2_size)[i]);
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner,
                                            &theta_grid_innercorner,
                                            (hydro_data->r0)[i]
                                                - 0.5*(hydro_data->r0_size)[i],
                                            (hydro_data->r1)[i]
                                                - 0.5*(hydro_data->r1_size)[i],
                                            0.0);
                hydroCoordinateToSpherical(&r_grid_outercorner,
                                            &theta_grid_outercorner,
                                            (hydro_data->r0)[i]
                                                + 0.5*(hydro_data->r0_size)[i],
                                            (hydro_data->r1)[i]
                                                + 0.5*(hydro_data->r1_size)[i],
                                            0.0);
            #endif
            int in_shell = (rmin <= r_grid_outercorner)
                        && (r_grid_innercorner  <  rmax)
                        && (theta_grid_outercorner >= theta_min)
                        && (theta_grid_innercorner  <  theta_max);

            if (in_shell)
            {
                double B = getMagneticFieldMagnitude(hydro_data, i);
                double V = hydroElementVolume(hydro_data, i);

                cell_index[j] = i;
                W_cell[j]     = (hydro_data->nonthermal_dens)[i] * B * B * V;
                j++;
            }
        }
    }

    /* ── Step 5: Weight-tuning loop  [mirrors photonEmitCyclosynch] ──────────
     *
     * The expected packet count for active cell j is:
     *   lambda_j = W_cell[j] / ph_weight_adjusted / fps
     *
     * analogous to:
     *   ph_dens_calc = integral_BB * V / ph_weight_adjusted / fps
     * in photonEmitCyclosynch.
     *
     * Iterate ph_weight_adjusted until:
     *   min_photons <= sum_j lambda_j <= max_photons
     *
     * Adjustment rules (identical to photonEmitCyclosynch):
     *   sum > max_photons : ph_weight_adjusted *= 10  (weight each packet
     *                        more, so fewer packets are needed)
     *   sum < min_photons : ph_weight_adjusted *= 0.5 (weight each packet
     *                        less, so more packets are needed)
     *
     * The loop operates on expected counts (before the Poisson draw) so
     * that the SYNCH_N_REF reference draws in buildSynchCellStrata are
     * not repeated during tuning.
     */
    double ph_weight_adjusted = ph_weight;
    double lambda_total       = 0.0;

    do
    {
        lambda_total = 0.0;
        for (j = 0; j < block_cnt; j++)
            lambda_total += W_cell[j] / ph_weight_adjusted / hydro_data->fps;

        if (lambda_total > (double)max_photons)
            ph_weight_adjusted *= 10.0;
        else if (lambda_total < (double)min_photons)
            ph_weight_adjusted *= 0.5;

    } while ((lambda_total > (double)max_photons) ||
             (lambda_total < (double)min_photons));

    fprintf(fPtr,
            ">> [photonEmitSynch] Converged: ph_weight_adjusted = %.3e, "
            "expected total packets = %.1f\n",
            ph_weight_adjusted, lambda_total);
    fflush(fPtr);

    /* ── Step 6: Poisson draw of per-cell packet counts ──────────────────────
     *
     * With ph_weight_adjusted converged, draw the actual packet count for
     * each active cell j from a Poisson distribution:
     *   ph_dens[j] ~ Poisson( W_cell[j] / ph_weight_adjusted / fps )
     *
     * The Poisson draw gives shot-noise-correct fluctuations around the
     * physical mean packet count, exactly as in photonEmitCyclosynch.
     * ph_tot accumulates the total across all cells.
     */
    int ph_tot = 0;
    for (j = 0; j < block_cnt; j++)
    {
        double lambda_j = W_cell[j] / ph_weight_adjusted / hydro_data->fps;
        ph_dens[j] = (int)gsl_ran_poisson(rand, lambda_j);
        ph_tot    += ph_dens[j];
    }

    fprintf(fPtr,
            ">> [photonEmitSynch] Poisson draw: ph_tot = %d packets "
            "across %d active cells.\n",
            ph_tot, block_cnt);
    fflush(fPtr);

    /* ── Step 7: Allocate and initialise output photon array ─────────────────
     *
     * Allocate ph_tot photon structs and initialise each to the canonical
     * null-photon state via initializePhoton. We fill these in the per-cell
     * loop in Step 8 and add them to the photon_list in a single batch call
     * in Step 9, mirroring the ph_emit pattern in photonEmitCyclosynch.
     */
    if (ph_tot == 0)
    {
        /*
         * All cells drew zero packets from the Poisson distribution.
         * This can happen when lambda_total is small and the Poisson
         * fluctuation rounds everything to zero. Log and return cleanly.
         */
        fprintf(fPtr,
                ">> [photonEmitSynch] WARNING: all cells drew ph_dens = 0 "
                "from the Poisson distribution (lambda_total = %.3f). "
                "No photons emitted.\n\n",
                lambda_total);
        fflush(fPtr);

        free(ph_dens);
        free(cell_index);
        free(W_cell);
        return 0;
    }

    struct photon *ph_emit = (struct photon *)
                             malloc(ph_tot * sizeof(struct photon));

    if (ph_emit == NULL)
    {
        fprintf(fPtr,
                ">> [photonEmitSynch] ERROR: malloc failed for ph_emit "
                "(%d photons). No photons emitted.\n\n", ph_tot);
        fflush(fPtr);

        free(ph_dens);
        free(cell_index);
        free(W_cell);
        return 0;
    }

    for (i = 0; i < ph_tot; i++)
        initializePhoton(&ph_emit[i]);

    /* ── Step 8: Per-cell emission loop ──────────────────────────────────────
     *
     * For each active cell j, if ph_dens[j] > 0:
     *
     *   (8a) Build per-cell strata from B_cell using buildSynchCellStrata.
     *        This constructs SYNCH_N_STRATA frequency intervals and their
     *        natural probabilities p_k anchored to B_cell, not a global
     *        representative field. Cells with very different B values
     *        (e.g. jet spine vs. sheath) therefore each get strata correctly
     *        centred on their own nu_c(B_cell).
     *
     *   (8b) Call emitCellPackets to fill ph_dens[j] photon structs into
     *        ph_emit[idx ... idx + ph_dens[j] - 1], advancing idx by
     *        ph_dens[j] on return.
     *
     *   (8c) Free the per-cell SynchCellStrata memory immediately to keep
     *        peak heap usage proportional to one cell's strata at a time.
     *
     * The previous implementation used `continue` to skip cells with
     * ph_dens[j] == 0. This is replaced by wrapping the cell body in
     * `if (ph_dens[j] > 0)` so cells with zero packets simply produce no
     * work without a jump statement.
     */
    int idx = 0;

    for (j = 0; j < block_cnt; j++)
    {
        if (ph_dens[j] > 0)
        {
            int ci = cell_index[j];

            /* ── (8a) Build per-cell strata anchored to B_cell ───────────── */
            SynchCellStrata cs;
            buildSynchCellStrata(&cs,
                                  getMagneticFieldMagnitude(hydro_data, ci),
                                  rand, fPtr);

            /* ── (8b) Emit ph_dens[j] packets from this cell ─────────────── */
            emitCellPackets(ci, ph_dens[j], &cs, ph_weight_adjusted,
                             ph_emit, &idx, hydro_data, rand, fPtr);

            /* ── (8c) Free per-cell strata ───────────────────────────────── */
            freeSynchCellStrata(&cs);
        }
    }

    /*
     * Sanity check: idx must equal ph_tot after the loop. A mismatch
     * indicates a logic error in the stratum allocation arithmetic inside
     * emitCellPackets. Log the discrepancy and truncate ph_tot to idx so
     * that addToPhotonList does not access uninitialised memory.
     */
    if (idx != ph_tot)
    {
        fprintf(fPtr,
                ">> [photonEmitSynch] WARNING: expected to fill %d photon "
                "structs but actually filled %d. "
                "Check stratum allocation logic in emitCellPackets.\n",
                ph_tot, idx);
        fflush(fPtr);
        ph_tot = idx;
    }

    /* ── Step 9: Add emitted photons to the photon list ──────────────────────
     *
     * A single batch call to addToPhotonList is more efficient than one
     * call per photon because it amortises the null-photon index scan over
     * the entire batch. addToPhotonList handles photon_list reallocation
     * if required and updates num_photons and num_null_photons correctly.
     * The ph_emit array is copied into photon_list->photons via memcpy
     * inside addToPhotonList, so ph_emit can be freed immediately after.
     */
    if (ph_tot > 0)
    {
        addToPhotonList(photon_list, ph_emit, ph_tot);
        n_emitted = ph_tot;
    }

    /* ── Step 10: Free all heap allocations ──────────────────────────────────
     *
     * ph_emit contents have been copied into photon_list by addToPhotonList
     * so freeing here is safe.
     */
    free(ph_emit);
    free(ph_dens);
    free(cell_index);
    free(W_cell);

    fprintf(fPtr,
            ">> [photonEmitSynch] Complete: emitted %d synchrotron photon "
            "packets from %d active cells "
            "(ph_weight_adjusted = %.3e).\n\n",
            n_emitted, block_cnt, ph_weight_adjusted);
    fflush(fPtr);

    return n_emitted;
}
