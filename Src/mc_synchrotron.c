/*
 * mc_synchrotron.c
 * ================
 * Pitch-angle averaged synchrotron photon emission with SSA weight
 * modification, integrated into the MCRaT photonList framework.
 *
 * References:
 *   Crusius & Schlickeiser (1986)  A&A 164, L16        [R(x) function]
 *   Ghisellini & Svensson (1991)   MNRAS 252, 313      [SSA kappa_nu]
 *   Kawashima et al. (2023)        ApJ 949, 101        [weight Eq. 40]
 */

#include "mc_synchrotron.h"

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 1: BESSEL INTEGRALS AND SPECTRAL FUNCTIONS
 * ══════════════════════════════════════════════════════════════════════════ */

static double bessel_K53_integrand(double xi, void *params)
{
    (void)params;
    if (xi > 500.0) return 0.0;
    return gsl_sf_bessel_Knu(5.0/3.0, xi);
}

/*
 * synchComputeRatX
 * ----------------
 * R(x) = F(x) = x * integral_x^inf K_{5/3}(xi) dxi
 * Uses GSL QAGIU on the semi-infinite interval.
 */
static double synchComputeRatX(double x, gsl_integration_workspace *ws)
{
    double result = 0.0, abserr = 0.0;
    gsl_function F;
    F.function = bessel_K53_integrand;
    F.params   = NULL;

    if (x > 50.0) return 0.0;

    int status = gsl_integration_qagiu(&F, x, 1e-14, 1e-10,
                                        1000, ws, &result, &abserr);
    if (status && x < 40.0)
    {
        /* non-fatal: QAGIU occasionally warns near machine precision */
    }
    return x * result;
}

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 2: UNIVERSAL TABLES
 * ══════════════════════════════════════════════════════════════════════════ */

void initSynchTables(SynchUniversalTables *tables, FILE *fPtr)
{
    int i, n_mono;
    double log_x_min = -5.0, log_x_max = 2.0;

    fprintf(fPtr, ">> [initSynchTables] Computing R(x) via Bessel integrals "
            "(%d points)...\n", SYNCH_N_X);
    fflush(fPtr);

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);

    /* ── R(x) array ──────────────────────────────────────────────────── */
    for (i = 0; i < SYNCH_N_X; i++)
    {
        double t        = (double)i / (SYNCH_N_X - 1);
        tables->x_arr[i] = pow(10.0, log_x_min + t*(log_x_max - log_x_min));
        tables->R_arr[i] = synchComputeRatX(tables->x_arr[i], ws);
    }
    gsl_integration_workspace_free(ws);

    /* ── Absorption kernel 2R + 2x dR/dx (central differences in log x) */
    for (i = 0; i < SYNCH_N_X; i++)
    {
        double dR_dlnx;
        if (i == 0)
            dR_dlnx = (tables->R_arr[1] - tables->R_arr[0]) /
                      (log(tables->x_arr[1]) - log(tables->x_arr[0]));
        else if (i == SYNCH_N_X - 1)
            dR_dlnx = (tables->R_arr[i] - tables->R_arr[i-1]) /
                      (log(tables->x_arr[i]) - log(tables->x_arr[i-1]));
        else
            dR_dlnx = (tables->R_arr[i+1] - tables->R_arr[i-1]) /
                      (log(tables->x_arr[i+1]) - log(tables->x_arr[i-1]));

        double ak = 2.0*tables->R_arr[i]
                  + 2.0*tables->x_arr[i] * (dR_dlnx / tables->x_arr[i]);
        tables->abs_kern_arr[i] = (ak > 0.0) ? ak : 0.0;
    }

    /* ── R(x) forward spline ─────────────────────────────────────────── */
    tables->R_acc    = gsl_interp_accel_alloc();
    tables->R_spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    gsl_spline_init(tables->R_spline,
                    tables->x_arr, tables->R_arr, SYNCH_N_X);

    /* ── abs_kern(x) forward spline ──────────────────────────────────── */
    tables->abs_kern_acc    = gsl_interp_accel_alloc();
    tables->abs_kern_spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    gsl_spline_init(tables->abs_kern_spline,
                    tables->x_arr, tables->abs_kern_arr, SYNCH_N_X);

    /* ── CDF of x ~ R(x) d(log x) ───────────────────────────────────── */
    tables->cdf_x[0] = 0.0;
    for (i = 1; i < SYNCH_N_X; i++)
    {
        double dlnx = log(tables->x_arr[i]) - log(tables->x_arr[i-1]);
        double f0   = tables->R_arr[i-1] * tables->x_arr[i-1];
        double f1   = tables->R_arr[i]   * tables->x_arr[i];
        tables->cdf_x[i] = tables->cdf_x[i-1] + 0.5*(f0 + f1)*dlnx;
    }
    double norm_x = tables->cdf_x[SYNCH_N_X - 1];
    for (i = 0; i < SYNCH_N_X; i++)
        tables->cdf_x[i] /= norm_x;

    /* Enforce strict monotonicity for inverse CDF spline */
    double cdf_mono[SYNCH_N_X], x_mono[SYNCH_N_X];
    n_mono = 0;
    for (i = 0; i < SYNCH_N_X; i++)
    {
        if (n_mono == 0 || tables->cdf_x[i] > cdf_mono[n_mono-1] + 1e-15)
        {
            cdf_mono[n_mono] = tables->cdf_x[i];
            x_mono[n_mono]   = tables->x_arr[i];
            n_mono++;
        }
    }
    tables->inv_F_acc    = gsl_interp_accel_alloc();
    tables->inv_F_spline = gsl_spline_alloc(gsl_interp_linear, n_mono);
    gsl_spline_init(tables->inv_F_spline, cdf_mono, x_mono, n_mono);

    /* ── CDF of alpha ~ (2/pi) sin^2(alpha) on [0, pi] ──────────────── */
    for (i = 0; i < SYNCH_N_ALPHA; i++)
    {
        double a              = (double)i / (SYNCH_N_ALPHA - 1) * M_PI;
        tables->alpha_arr[i]  = a;
        tables->cdf_alpha[i]  = (a - 0.5*sin(2.0*a)) / M_PI;
    }
    double cdf_am[SYNCH_N_ALPHA], alpha_am[SYNCH_N_ALPHA];
    n_mono = 0;
    for (i = 0; i < SYNCH_N_ALPHA; i++)
    {
        if (n_mono == 0 ||
            tables->cdf_alpha[i] > cdf_am[n_mono-1] + 1e-15)
        {
            cdf_am[n_mono]   = tables->cdf_alpha[i];
            alpha_am[n_mono] = tables->alpha_arr[i];
            n_mono++;
        }
    }
    tables->inv_alpha_acc    = gsl_interp_accel_alloc();
    tables->inv_alpha_spline = gsl_spline_alloc(gsl_interp_linear, n_mono);
    gsl_spline_init(tables->inv_alpha_spline, cdf_am, alpha_am, n_mono);

    fprintf(fPtr, ">> [initSynchTables] All tables ready.\n");
    fflush(fPtr);
}

void freeSynchTables(SynchUniversalTables *tables)
{
    gsl_spline_free(tables->R_spline);
    gsl_spline_free(tables->abs_kern_spline);
    gsl_spline_free(tables->inv_F_spline);
    gsl_spline_free(tables->inv_alpha_spline);
    gsl_interp_accel_free(tables->R_acc);
    gsl_interp_accel_free(tables->abs_kern_acc);
    gsl_interp_accel_free(tables->inv_F_acc);
    gsl_interp_accel_free(tables->inv_alpha_acc);
}

/* ── Inline samplers ─────────────────────────────────────────────────── */

static inline double synchSampleX(const SynchUniversalTables *t, double u)
{
    return gsl_spline_eval(t->inv_F_spline, u, t->inv_F_acc);
}

static inline double synchSampleAlpha(const SynchUniversalTables *t, double u)
{
    return gsl_spline_eval(t->inv_alpha_spline, u, t->inv_alpha_acc);
}

static inline double synchEvalAbsKern(const SynchUniversalTables *t, double x)
{
    if (x <= t->x_arr[0])             return t->abs_kern_arr[0];
    if (x >= t->x_arr[SYNCH_N_X - 1]) return 0.0;
    return gsl_spline_eval(t->abs_kern_spline, x, t->abs_kern_acc);
}

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 3: GAMMA SAMPLER
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * synchSampleGammaEmission
 * ------------------------
 * Sample gamma ~ N(gamma)*gamma^2 ∝ gamma^{2-p} for synchrotron EMISSION.
 * This is the correct weighting for photon emission since each electron's
 * contribution to the emissivity scales as gamma^2.
 *
 * This function is kept as the bespoke analytic inverse CDF because the
 * existing electron.c samplers (samplePowerLaw, sampleBrokenPowerLawSubgroup)
 * sample from N(gamma) directly, which is the wrong weighting for emission.
 * They are used instead in synchEmitOneNuRef() for the reference CDF only.
 */
static double synchSampleGammaEmission(gsl_rng *rand)
{
    double u = gsl_rng_uniform_pos(rand);

    #if NONTHERMAL_E_DIST == POWERLAW

        double p    = POWERLAW_INDEX;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;
        double q    = 3.0 - p;    /* exponent for gamma^{2-p} weighting */

        if (fabs(q) < 1e-6)
            return gmin * pow(gmax/gmin, u);   /* p==3: log-uniform in gamma */

        double glo_q = pow(gmin, q);
        double ghi_q = pow(gmax, q);
        return pow(glo_q + u*(ghi_q - glo_q), 1.0/q);

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        double p1   = POWERLAW_INDEX_1;
        double p2   = POWERLAW_INDEX_2;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;
        double gbr  = GAMMA_BREAK;
        double q1   = 3.0 - p1;
        double q2   = 3.0 - p2;

        /* Continuity factor at the break for the gamma^2-weighted distribution */
        double C_cont = pow(gbr, p2 - p1);

        /* Power (integral of N(gamma)*gamma^2) in each segment */
        double I1 = (fabs(q1) < 1e-6)
                    ? log(gbr/gmin)
                    : (pow(gbr,q1) - pow(gmin,q1)) / q1;
        double I2 = C_cont * ((fabs(q2) < 1e-6)
                               ? log(gmax/gbr)
                               : (pow(gmax,q2) - pow(gbr,q2)) / q2);
        double f1 = I1 / (I1 + I2);

        if (u < f1)
        {
            double u1 = u / f1;
            if (fabs(q1) < 1e-6) return gmin * pow(gbr/gmin, u1);
            double glo_q = pow(gmin, q1);
            double ghi_q = pow(gbr,  q1);
            return pow(glo_q + u1*(ghi_q - glo_q), 1.0/q1);
        }
        else
        {
            double u2 = (u - f1) / (1.0 - f1);
            if (fabs(q2) < 1e-6) return gbr * pow(gmax/gbr, u2);
            double glo_q = pow(gbr,  q2);
            double ghi_q = pow(gmax, q2);
            return pow(glo_q + u2*(ghi_q - glo_q), 1.0/q2);
        }

    #endif
}

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 4: KAPPA_NU TABLE
 * ══════════════════════════════════════════════════════════════════════════ */


SynchKappaTable *buildSynchKappaTable(double B,
                                       const SynchUniversalTables *tables,
                                       FILE *fPtr)
{
    /*
     * Build the opacity table with n_e = 1 (unit number density).
     * The actual cell n_e is applied at lookup time in synchKappaAtNuScaled.
     *
     * N(gamma) is evaluated using singleElectronPowerLaw /
     * singleElectronBrokenPowerLaw from electron.c, which return the
     * normalised shape hat_N(gamma) = N(gamma)/n_e.  This keeps all
     * normalisation logic in one place and avoids duplicating it here.
     */
    int j, k;
    SynchKappaTable *kt = (SynchKappaTable *)malloc(sizeof(SynchKappaTable));
    kt->n_nu      = SYNCH_N_NU;
    kt->log_nu    = (double *)malloc(SYNCH_N_NU * sizeof(double));
    kt->log_kappa = (double *)malloc(SYNCH_N_NU * sizeof(double));
    kt->B_ref     = B;
    kt->ne_ref    = 1.0;   /* table is built for n_e = 1; scale at lookup */

    double nu_s_coeff = (3.0*CHARGE_EL) / (4.0*M_PI*M_EL*C_LIGHT);
    double nu_lo = 1e-4 * nu_s_coeff * B * GAMMA_MIN*GAMMA_MIN;
    double nu_hi = 1e4  * nu_s_coeff * B * GAMMA_MAX*GAMMA_MAX;

    /* Prefactor: sqrt(3) e^3 B / (8 pi me c)  [Ghisellini & Svensson 1991] */
    double prefac = (sqrt(3.0)*CHARGE_EL*CHARGE_EL*CHARGE_EL * B)
                  / (8.0*M_PI*M_EL*C_LIGHT);

    /* Gauss-Legendre nodes in log(gamma) */
    gsl_integration_glfixed_table *gl =
        gsl_integration_glfixed_table_alloc(SYNCH_N_GL);

    for (j = 0; j < SYNCH_N_NU; j++)
    {
        double t  = (double)j / (SYNCH_N_NU - 1);
        double nu = nu_lo * pow(nu_hi/nu_lo, t);
        kt->log_nu[j] = log10(nu);

        double integral = 0.0;
        for (k = 0; k < SYNCH_N_GL; k++)
        {
            double xi, wi;
            gsl_integration_glfixed_point(log(GAMMA_MIN), log(GAMMA_MAX),
                                           k, &xi, &wi, gl);
            double gam = exp(xi);
            double x_k = nu / (nu_s_coeff * B * gam*gam);
            double ak  = synchEvalAbsKern(tables, x_k);

            /*
             * hat_N(gamma) = N(gamma)/n_e evaluated using electron.c:
             *   singleElectronPowerLaw    -> A * gamma^{-p}
             *   singleElectronBrokenPowerLaw -> A * gamma^{-p1 or -p2}
             * Both functions already include the normalisation constant A.
             * Multiply by gamma (Jacobian for d log gamma measure).
             */
            #if NONTHERMAL_E_DIST == POWERLAW
                double Nhat_x_gam = singleElectronPowerLaw(gam,
                                                            POWERLAW_INDEX,
                                                            GAMMA_MIN,
                                                            GAMMA_MAX) * gam;
            #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
                double Nhat_x_gam = singleElectronBrokenPowerLaw(gam,
                                                                   POWERLAW_INDEX_1,
                                                                   POWERLAW_INDEX_2,
                                                                   GAMMA_MIN,
                                                                   GAMMA_MAX,
                                                                   GAMMA_BREAK) * gam;
            #endif
            /* n_e = 1 here; actual n_e applied at lookup via scaling */
            integral += wi * Nhat_x_gam * ak;
        }

        double kap = (prefac / (nu*nu)) * integral;
        kt->log_kappa[j] = log10((kap > 1e-300) ? kap : 1e-300);
    }

    gsl_integration_glfixed_table_free(gl);

    kt->acc    = gsl_interp_accel_alloc();
    kt->spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_NU);
    gsl_spline_init(kt->spline, kt->log_nu, kt->log_kappa, SYNCH_N_NU);

    fprintf(fPtr,
            ">> [buildSynchKappaTable] Built unit kappa_nu table: "
            "B=%.2e G  ne=1  nu=[%.2e, %.2e] Hz\n",
            B, pow(10.0, kt->log_nu[0]),
            pow(10.0, kt->log_nu[SYNCH_N_NU-1]));
    fflush(fPtr);

    return kt;
}

void freeSynchKappaTable(SynchKappaTable *kt)
{
    gsl_spline_free(kt->spline);
    gsl_interp_accel_free(kt->acc);
    free(kt->log_nu);
    free(kt->log_kappa);
    free(kt);
}

/*
 * synchKappaAtNu
 * --------------
 * Evaluate kappa_nu [cm^{-1}] at frequency nu [Hz] via log-log interpolation.
 * The table was built for (B_ref, ne_ref); rescale linearly for other cells:
 *   kappa_nu(B, ne) = kappa_nu_table * (ne/ne_ref) * (B/B_ref)
 * The B scaling follows from the prefactor sqrt(3) e^3 B / (8 pi me c).
 * Pass B_cell=kt->B_ref and ne_cell=kt->ne_ref to use the table as-is.
 */
static double synchKappaAtNu(double nu, const SynchKappaTable *kt)
{
    double lnu = log10(nu);
    if (lnu <= kt->log_nu[0])
        return pow(10.0, kt->log_kappa[0]);
    if (lnu >= kt->log_nu[kt->n_nu - 1])
        return pow(10.0, kt->log_kappa[kt->n_nu - 1]);
    return pow(10.0, gsl_spline_eval(kt->spline, lnu, kt->acc));
}

/*
 * synchKappaAtNuScaled
 * --------------------
 * Evaluate the physical SSA opacity kappa_nu [cm^{-1}] for a grid cell
 * with magnetic field B_cell [G] and nonthermal electron number density
 * ne_cell [cm^{-3}].
 *
 * The kappa table was built with n_e = 1 and B = B_ref, so the full
 * physical opacity is recovered by rescaling linearly:
 *
 *   kappa_nu(B_cell, ne_cell) = kappa_table(nu; B_ref, ne=1)
 *                               * ne_cell          [n_e scaling]
 *                               * (B_cell / B_ref) [B scaling]
 *
 * Both scalings are exact consequences of the linear dependence of
 * kappa_nu on n_e and B in the Ghisellini & Svensson (1991) expression.
 *
 * This is the ONLY function that should be called externally to obtain
 * a physical opacity. synchKappaAtNu is now private (static).
 */
static inline double synchKappaAtNuScaled(double nu,
                                           const SynchKappaTable *kt,
                                           double B_cell,
                                           double ne_cell)
{
    return synchKappaAtNu(nu, kt)
           * ne_cell
           * (B_cell / kt->B_ref);
}

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 5: STRATIFIED FREQUENCY SAMPLER SETUP
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * synchEmitOneNu
 * --------------
 * Draw a single representative photon frequency from the natural emission
 * distribution of cell i. Used ONLY to build the reference CDF in
 * buildSynchStratifiedParams — not for actual photon emission.
 *
 * Here we use samplePowerLaw / sampleBrokenPowerLawSubgroup from electron.c
 * to sample gamma ~ N(gamma). This slightly mis-weights the reference
 * distribution relative to the true N(gamma)*gamma^2 emission weighting,
 * but since the reference CDF only needs to identify which frequency strata
 * are populated (not their exact relative weights), this is sufficient and
 * avoids duplicating sampling logic.
 */
static double synchEmitOneNu(int i,
                               struct hydro_dataframe *hydro_data,
                               const SynchUniversalTables *tables,
                               double b_field,
                               gsl_rng *rand)
{
    double nu_s_coeff = (3.0*CHARGE_EL) / (4.0*M_PI*M_EL*C_LIGHT);

    /* Sample gamma from N(gamma) using existing electron.c functions.
     * These are not weighted by gamma^2 but are correct for identifying
     * which frequency strata receive emission. */
    #if NONTHERMAL_E_DIST == POWERLAW
        double gamma_k = samplePowerLaw(POWERLAW_INDEX,
                                         GAMMA_MIN, GAMMA_MAX,
                                         rand, NULL);
    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
        double gamma_k = sampleBrokenPowerLawSubgroup(POWERLAW_INDEX_1,
                                                       POWERLAW_INDEX_2,
                                                       GAMMA_MIN, GAMMA_MAX,
                                                       GAMMA_BREAK,
                                                       rand, NULL);
    #endif

    double alpha_k = synchSampleAlpha(tables, gsl_rng_uniform_pos(rand));
    double x_k     = synchSampleX(tables,     gsl_rng_uniform_pos(rand));
    double nu_c    = nu_s_coeff * b_field * sin(alpha_k) * gamma_k*gamma_k;
    return x_k * nu_c;
}

void buildSynchStratifiedParams(SynchStratifiedParams *sp,
                                 const SynchUniversalTables *tables,
                                 struct hydro_dataframe *hydro_data,
                                 double nu_min, double nu_max,
                                 gsl_rng *rand, FILE *fPtr)
{
    int i, k, n_mono, bin;
    double log_nu_min = log10(nu_min);
    double log_nu_max = log10(nu_max);
    double dlog       = (log_nu_max - log_nu_min) / SYNCH_N_CDF_BINS;

    sp->n_strata = SYNCH_N_STRATA;

    /* Log-spaced stratum edges */
    for (k = 0; k <= SYNCH_N_STRATA; k++)
    {
        double t = (double)k / SYNCH_N_STRATA;
        sp->strata_edges[k] = pow(10.0,
                                   log_nu_min + t*(log_nu_max - log_nu_min));
    }

    /* ── Draw reference sample from all emitting cells ───────────────── */
    fprintf(fPtr,
            ">> [buildSynchStratifiedParams] Drawing %d reference photons "
            "to measure stratum probabilities...\n", SYNCH_N_REF);
    fflush(fPtr);

    /* Weight cells by B^2 * ne to select emitting cell for each ref photon */
    double *cell_weights = (double *)malloc(hydro_data->num_elements
                                             * sizeof(double));
    double  wsum = 0.0;
    for (i = 0; i < hydro_data->num_elements; i++)
    {
        double b = getMagneticFieldMagnitude(hydro_data, i);
        /* nonthermal_dens is defined when NONTHERMAL_E_DIST != OFF */
        double ne_i = (hydro_data->nonthermal_dens)[i];
        cell_weights[i] = b*b * ne_i;
        wsum += cell_weights[i];
    }
    /* Build a simple alias table for cell selection */
    /* Re-use the alias logic from the existing codebase style:
       walk the CDF with a binary search for simplicity here since
       this is only called once at initialisation                   */
    double *cell_cdf = (double *)malloc(hydro_data->num_elements
                                         * sizeof(double));
    cell_cdf[0] = cell_weights[0] / wsum;
    for (i = 1; i < hydro_data->num_elements; i++)
        cell_cdf[i] = cell_cdf[i-1] + cell_weights[i] / wsum;

    /* Histogram counters */
    int stratum_counts[SYNCH_N_STRATA];
    memset(stratum_counts, 0, sizeof(stratum_counts));
    int hist[SYNCH_N_CDF_BINS];
    memset(hist, 0, sizeof(hist));

    for (i = 0; i < SYNCH_N_REF; i++)
    {
        /* Select a cell proportional to B^2 * ne */
        double u_cell = gsl_rng_uniform_pos(rand);
        int    cell   = 0;
        /* Binary search for cell index in CDF */
        int lo = 0, hi = hydro_data->num_elements - 1;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (cell_cdf[mid] < u_cell) lo = mid + 1;
            else                         hi = mid;
        }
        cell = lo;

        double b = getMagneticFieldMagnitude(hydro_data, cell);
        double nu = synchEmitOneNu(cell, hydro_data, tables, b, rand);

        /* Clamp to global range */
        if (nu < nu_min) nu = nu_min;
        if (nu > nu_max) nu = nu_max * (1.0 - 1e-10);

        /* Stratum count */
        for (k = 0; k < SYNCH_N_STRATA; k++)
        {
            if (nu >= sp->strata_edges[k] && nu < sp->strata_edges[k+1])
            {
                stratum_counts[k]++;
                break;
            }
        }

        /* Fine CDF histogram */
        bin = (int)((log10(nu) - log_nu_min) / dlog);
        if (bin >= 0 && bin < SYNCH_N_CDF_BINS) hist[bin]++;
    }

    free(cell_weights);
    free(cell_cdf);

    /* Stratum probabilities */
    fprintf(fPtr, ">> [buildSynchStratifiedParams] Stratum probabilities:\n");
    for (k = 0; k < SYNCH_N_STRATA; k++)
    {
        sp->stratum_probs[k] = (double)stratum_counts[k] / SYNCH_N_REF;
        fprintf(fPtr, "     stratum %2d  [%.2e, %.2e] Hz  p = %.6f\n",
                k, sp->strata_edges[k], sp->strata_edges[k+1],
                sp->stratum_probs[k]);
    }
    fflush(fPtr);

    /* ── Build marginal CDF of log10(nu) ─────────────────────────────── */
    double cum = 0.0;
    sp->cdf_log_nu_edges[0] = log_nu_min;
    sp->cdf_log_nu_vals[0]  = 0.0;
    for (i = 0; i < SYNCH_N_CDF_BINS; i++)
    {
        cum += (double)hist[i];
        sp->cdf_log_nu_edges[i+1] = log_nu_min + (i+1)*dlog;
        sp->cdf_log_nu_vals[i+1]  = cum / SYNCH_N_REF;
    }
    sp->cdf_log_nu_vals[SYNCH_N_CDF_BINS] = 1.0; /* enforce endpoint */

    /* Enforce strict monotonicity */
    double mono_cdf[SYNCH_N_CDF_BINS+1], mono_edges[SYNCH_N_CDF_BINS+1];
    n_mono = 0;
    for (i = 0; i <= SYNCH_N_CDF_BINS; i++)
    {
        if (n_mono == 0 ||
            sp->cdf_log_nu_vals[i] > mono_cdf[n_mono-1] + 1e-12)
        {
            mono_cdf[n_mono]   = sp->cdf_log_nu_vals[i];
            mono_edges[n_mono] = sp->cdf_log_nu_edges[i];
            n_mono++;
        }
    }
    sp->inv_nu_cdf_acc    = gsl_interp_accel_alloc();
    sp->inv_nu_cdf_spline = gsl_spline_alloc(gsl_interp_linear, n_mono);
    gsl_spline_init(sp->inv_nu_cdf_spline, mono_cdf, mono_edges, n_mono);

    fprintf(fPtr,
            ">> [buildSynchStratifiedParams] Stratified sampler ready.\n");
    fflush(fPtr);
}

void freeSynchStratifiedParams(SynchStratifiedParams *sp)
{
    gsl_spline_free(sp->inv_nu_cdf_spline);
    gsl_interp_accel_free(sp->inv_nu_cdf_acc);
}

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 6: SINGLE-PHOTON EMISSION HELPER
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * synchFillPhoton
 * ---------------
 * Populate all fields of a single struct photon for synchrotron emission
 * from hydro cell i at comoving frequency fr_dum.
 *
 * Steps:
 *  1. Build isotropic comoving 4-momentum at fr_dum
 *  2. Lorentz-boost to lab frame (mirrors photonEmitCyclosynch)
 *  3. Assign random position within cell
 *  4. Set Stokes, weight, type, and bookkeeping fields
 */
static void synchFillPhoton(struct photon *ph,
                              int            i,
                              double         fr_dum,
                              double         ph_weight_adjusted,
                              struct hydro_dataframe *hydro_data,
                              gsl_rng       *rand,
                              FILE          *fPtr)
{
    double com_v_phi, com_v_theta;
    double position_phi = 0.0;
    double position_rand = 0.0, position2_rand = 0.0, position3_rand = 0.0;
    double cartesian_position_rand_array[3];
    double *p_comv = NULL, *boost = NULL, *l_boost = NULL;

    p_comv  = (double *)malloc(4 * sizeof(double));
    boost   = (double *)malloc(3 * sizeof(double));
    l_boost = (double *)malloc(4 * sizeof(double));

    /* ── (1) Random isotropic direction in comoving frame ────────────── */
    /*     Matches the same convention as photonEmitCyclosynch           */
    com_v_phi   = gsl_rng_uniform(rand) * 2.0 * M_PI;
    com_v_theta = acos((gsl_rng_uniform(rand) * 2.0) - 1.0);

    *(p_comv+0) = PL_CONST * fr_dum / C_LIGHT;
    *(p_comv+1) = (PL_CONST * fr_dum / C_LIGHT) * sin(com_v_theta) * cos(com_v_phi);
    *(p_comv+2) = (PL_CONST * fr_dum / C_LIGHT) * sin(com_v_theta) * sin(com_v_phi);
    *(p_comv+3) = (PL_CONST * fr_dum / C_LIGHT) * cos(com_v_theta);

    /* ── (2) Lorentz boost to lab frame ──────────────────────────────── */
    /*     Mirrors the convention in photonEmitCyclosynch exactly        */
#if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
    position_phi = gsl_rng_uniform(rand) * 2.0 * M_PI;
#else
    position_phi = 0.0;
#endif

#if DIMENSIONS == THREE
    hydroVectorToCartesian(boost,
                            (hydro_data->v0)[i],
                            (hydro_data->v1)[i],
                            (hydro_data->v2)[i],
                            (hydro_data->r0)[i],
                            (hydro_data->r1)[i],
                            (hydro_data->r2)[i]);
#elif DIMENSIONS == TWO_POINT_FIVE
    hydroVectorToCartesian(boost,
                            (hydro_data->v0)[i],
                            (hydro_data->v1)[i],
                            (hydro_data->v2)[i],
                            (hydro_data->r0)[i],
                            (hydro_data->r1)[i],
                            position_phi);
#else
    hydroVectorToCartesian(boost,
                            (hydro_data->v0)[i],
                            (hydro_data->v1)[i],
                            0,
                            (hydro_data->r0)[i],
                            (hydro_data->r1)[i],
                            position_phi);
#endif
    /* Negate boost to go from fluid -> lab frame (same sign convention  */
    /* as photonEmitCyclosynch)                                          */
    (*(boost+0)) *= -1.0;
    (*(boost+1)) *= -1.0;
    (*(boost+2)) *= -1.0;

    lorentzBoost(boost, p_comv, l_boost, 'p', fPtr);

    ph->p0     = *(l_boost+0);
    ph->p1     = *(l_boost+1);
    ph->p2     = *(l_boost+2);
    ph->p3     = *(l_boost+3);
    ph->comv_p0 = *(p_comv+0);
    ph->comv_p1 = *(p_comv+1);
    ph->comv_p2 = *(p_comv+2);
    ph->comv_p3 = *(p_comv+3);

    /* ── (3) Random position within cell ─────────────────────────────── */
    position_rand  = gsl_rng_uniform_pos(rand) * (hydro_data->r0_size)[i]
                   - 0.5 * (hydro_data->r0_size)[i];
    position2_rand = gsl_rng_uniform_pos(rand) * (hydro_data->r1_size)[i]
                   - 0.5 * (hydro_data->r1_size)[i];

#if DIMENSIONS == THREE
    position3_rand = gsl_rng_uniform_pos(rand) * (hydro_data->r2_size)[i]
                   - 0.5 * (hydro_data->r2_size)[i];
    hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array,
                                      (hydro_data->r0)[i] + position_rand,
                                      (hydro_data->r1)[i] + position2_rand,
                                      (hydro_data->r2)[i] + position3_rand);
#else
    hydroCoordinateToMcratCoordinate(&cartesian_position_rand_array,
                                      (hydro_data->r0)[i] + position_rand,
                                      (hydro_data->r1)[i] + position2_rand,
                                      position_phi);
#endif

    ph->r0 = cartesian_position_rand_array[0];
    ph->r1 = cartesian_position_rand_array[1];
    ph->r2 = cartesian_position_rand_array[2];

    /* ── (4) Stokes, weight, bookkeeping ─────────────────────────────── */
    ph->s0  = 1;  /* unpolarised: I=1, Q=U=V=0 */
    ph->s1  = 0;
    ph->s2  = 0;
    ph->s3  = 0;
    ph->num_scatt          = 0;
    ph->weight             = ph_weight_adjusted;
    ph->nearest_block_index = i;
    ph->type               = SYNCH_PHOTON;
    ph->recalc_properties  = 1;

    free(p_comv); free(boost); free(l_boost);
}

/* ══════════════════════════════════════════════════════════════════════════
 * SECTION 7: MAIN EMISSION FUNCTION
 * ══════════════════════════════════════════════════════════════════════════ */

int photonEmitSynch(struct photonList          *photon_list,
                    double                      r_inj,
                    double                      ph_weight,
                    int                         maximum_photons,
                    double                      theta_min,
                    double                      theta_max,
                    struct hydro_dataframe     *hydro_data,
                    const SynchUniversalTables *tables,
                    const SynchStratifiedParams *sp,
                    gsl_rng                    *rand,
                    FILE                       *fPtr)
{
    int    i, k, j;
    int    block_cnt = 0, ph_tot = 0;
    double r_grid_innercorner = 0, r_grid_outercorner   = 0;
    double theta_grid_innercorner = 0, theta_grid_outercorner = 0;
    double rmin = 0, rmax = 0;
    double b_field = 0, ne_cell = 0;
    double nu_s_coeff = (3.0*CHARGE_EL) / (4.0*M_PI*M_EL*C_LIGHT);

    /* ── Step 1: Determine injection radial limits ───────────────────── */
    /*   Mirrors photonEmitCyclosynch: use the light-travel window       */
    rmin = calcCyclosynchRLimits(hydro_data->scatt_frame_number,
                                  hydro_data->inj_frame_number,
                                  hydro_data->fps, r_inj, "min");
    rmax = calcCyclosynchRLimits(hydro_data->scatt_frame_number,
                                  hydro_data->inj_frame_number,
                                  hydro_data->fps, r_inj, "max");

    fprintf(fPtr,
            ">> [photonEmitSynch] rmin=%.3e rmax=%.3e "
            "theta=[%.3f, %.3f] rad\n",
            rmin, rmax, theta_min, theta_max);
    fflush(fPtr);

    /* ── Step 2: Count eligible cells and compute emission weights ───── */
    for (i = 0; i < hydro_data->num_elements; i++)
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
                                        0);
            hydroCoordinateToSpherical(&r_grid_outercorner,
                                        &theta_grid_outercorner,
                                        (hydro_data->r0)[i]
                                            + 0.5*(hydro_data->r0_size)[i],
                                        (hydro_data->r1)[i]
                                            + 0.5*(hydro_data->r1_size)[i],
                                        0);
        #endif
        if ((rmin       <= r_grid_outercorner)   &&
            (r_grid_innercorner  < rmax)          &&
            (theta_grid_outercorner >= theta_min) &&
            (theta_grid_innercorner  < theta_max))
        {
            block_cnt++;
        }
    }

    fprintf(fPtr,
            ">> [photonEmitSynch] %d hydro elements selected for emission.\n",
            block_cnt);
    fflush(fPtr);

    //TODO: make this throw an error
    if (block_cnt == 0) return 0;

    /* ── Step 3: Determine photons per stratum ───────────────────────── */
    /*   Cap total at maximum_photons; distribute evenly across strata   */
    int n_active_strata = 0;
    double P_MIN = 1e-7;
    for (k = 0; k < SYNCH_N_STRATA; k++)
        if (sp->stratum_probs[k] >= P_MIN) n_active_strata++;

    if (n_active_strata == 0)
    {
        fprintf(fPtr, ">> [photonEmitSynch] WARNING: no active strata.\n");
        fflush(fPtr);
        return 0;
    }

    int n_per_stratum = maximum_photons / n_active_strata;
    if (n_per_stratum < 1) n_per_stratum = 1;

    int n_total_emit = n_per_stratum * n_active_strata;

    /* Allocate temporary array for all new photons */
    struct photon *ph_emit = (struct photon *)malloc(
                                n_total_emit * sizeof(struct photon));
    if (!ph_emit)
    {
        fprintf(fPtr, ">> [photonEmitSynch] ERROR: malloc failed.\n");
        fflush(fPtr);
        return 0;
    }

    /* ── Step 4: Build cell selection alias (B^2 * ne weights) ──────── */
    /*   Only over eligible cells; reuse simple CDF here                 */
    int    *eligible    = (int    *)malloc(block_cnt * sizeof(int));
    double *cell_cdf    = (double *)malloc(block_cnt * sizeof(double));
    int     el_cnt      = 0;
    double  wsum        = 0.0;

    for (i = 0; i < hydro_data->num_elements; i++)
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
                                        0);
            hydroCoordinateToSpherical(&r_grid_outercorner,
                                        &theta_grid_outercorner,
                                        (hydro_data->r0)[i]
                                            + 0.5*(hydro_data->r0_size)[i],
                                        (hydro_data->r1)[i]
                                            + 0.5*(hydro_data->r1_size)[i],
                                        0);
        #endif
        if ((rmin       <= r_grid_outercorner)   &&
            (r_grid_innercorner  < rmax)          &&
            (theta_grid_outercorner >= theta_min) &&
            (theta_grid_innercorner  < theta_max))
        {
            eligible[el_cnt] = i;
            b_field  = getMagneticFieldMagnitude(hydro_data, i);
            ne_cell  = (hydro_data->nonthermal_dens)[i];
            double w = b_field*b_field * ne_cell
                     * hydroElementVolume(hydro_data, i);
            cell_cdf[el_cnt] = (el_cnt == 0) ? w : cell_cdf[el_cnt-1] + w;
            wsum += w;
            el_cnt++;
        }
    }
    /* Normalise CDF */
    for (i = 0; i < el_cnt; i++) cell_cdf[i] /= wsum;

    /* ── Step 5: Emit photons stratum by stratum ──────────────────────── */
    double uniform_prob = 1.0 / SYNCH_N_STRATA;
    ph_tot = 0;

    for (k = 0; k < SYNCH_N_STRATA; k++)
    {
        double p_k = sp->stratum_probs[k];

        if (p_k >= P_MIN)
        {
            double w_correction      = p_k / uniform_prob;
            double ph_weight_adjusted = ph_weight * w_correction;

            double log_nu_lo = log10(sp->strata_edges[k]);
            double log_nu_hi = log10(sp->strata_edges[k+1]);

            /* Forward CDF at lower stratum edge */
            double C_lo = 0.0;
            for (i = 1; i <= SYNCH_N_CDF_BINS; i++)
            {
                if (sp->cdf_log_nu_edges[i] >= log_nu_lo)
                {
                    double t = (log_nu_lo - sp->cdf_log_nu_edges[i-1])
                             / (sp->cdf_log_nu_edges[i]
                                - sp->cdf_log_nu_edges[i-1]);
                    C_lo = sp->cdf_log_nu_vals[i-1]
                         + t*(sp->cdf_log_nu_vals[i]
                              - sp->cdf_log_nu_vals[i-1]);
                    break;
                }
            }

            /* Forward CDF at upper stratum edge */
            double C_hi = 1.0;
            for (i = 1; i <= SYNCH_N_CDF_BINS; i++)
            {
                if (sp->cdf_log_nu_edges[i] >= log_nu_hi)
                {
                    double t = (log_nu_hi - sp->cdf_log_nu_edges[i-1])
                             / (sp->cdf_log_nu_edges[i]
                                - sp->cdf_log_nu_edges[i-1]);
                    C_hi = sp->cdf_log_nu_vals[i-1]
                         + t*(sp->cdf_log_nu_vals[i]
                              - sp->cdf_log_nu_vals[i-1]);
                    break;
                }
            }

            double dC = C_hi - C_lo;

            if (dC >= 1e-10)
            {
                fprintf(fPtr,
                        ">> [photonEmitSynch] stratum %2d  "
                        "[%.2e, %.2e] Hz  p=%.4f  w_corr=%.4f  N=%d\n",
                        k, sp->strata_edges[k], sp->strata_edges[k+1],
                        p_k, w_correction, n_per_stratum);
                fflush(fPtr);

                for (j = 0; j < n_per_stratum; j++)
                {
                    /* (a) Select emitting cell from eligible list via CDF */
                    double u_cell = gsl_rng_uniform_pos(rand);
                    int    lo     = 0, hi = el_cnt - 1;
                    while (lo < hi)
                    {
                        int mid = (lo + hi) / 2;
                        if (cell_cdf[mid] < u_cell) lo = mid + 1;
                        else                         hi = mid;
                    }
                    int cell_idx = eligible[lo];

                    /* (b) Get local B and nonthermal electron density */
                    b_field = getMagneticFieldMagnitude(hydro_data, cell_idx);
                    ne_cell = (hydro_data->nonthermal_dens)[cell_idx];

                    /* (c) Sample gamma and alpha from their natural
                     *     distributions */
                    double gamma_k = synchSampleGammaEmission(rand);
                    double alpha_k = synchSampleAlpha(tables,
                                                       gsl_rng_uniform_pos(rand));

                    /* (d) Conditional nu sampling via inverse marginal CDF */
                    double u_cond  = C_lo + gsl_rng_uniform_pos(rand) * dC;
                    double u_lo_sp = sp->cdf_log_nu_vals[0];
                    double u_hi_sp = sp->cdf_log_nu_vals[SYNCH_N_CDF_BINS];
                    if (u_cond < u_lo_sp + 1e-12) u_cond = u_lo_sp + 1e-12;
                    if (u_cond > u_hi_sp - 1e-12) u_cond = u_hi_sp - 1e-12;

                    double log_nu_k = gsl_spline_eval(sp->inv_nu_cdf_spline,
                                                       u_cond,
                                                       sp->inv_nu_cdf_acc);
                    double fr_dum = pow(10.0, log_nu_k);

                    /* Clamp strictly within stratum */
                    if (fr_dum < sp->strata_edges[k])
                        fr_dum = sp->strata_edges[k];
                    if (fr_dum > sp->strata_edges[k+1])
                        fr_dum = sp->strata_edges[k+1];

                    /* (e) Fill the photon struct */
                    synchFillPhoton(&ph_emit[ph_tot],
                                     cell_idx,
                                     fr_dum,
                                     ph_weight_adjusted,
                                     hydro_data,
                                     rand,
                                     fPtr);
                    ph_tot++;
                }
            }
            else
            {
                fprintf(fPtr,
                        ">> [photonEmitSynch] WARNING: stratum %d passed "
                        "P_MIN but has degenerate CDF interval (dC=%.2e). "
                        "Skipping — increase SYNCH_N_REF if this recurs.\n",
                        k, dC);
                fflush(fPtr);
            }
        }
    }

    free(eligible);
    free(cell_cdf);

    /* ── Step 6: Add all emitted photons to the photon list ───────────── */
    /*   addToPhotonList handles memory growth and null-slot reuse,        */
    /*   exactly as photonEmitCyclosynch does via setPhotonList            */
    addToPhotonList(photon_list, ph_emit, (size_t)ph_tot);
    free(ph_emit);

    fprintf(fPtr,
            ">> [photonEmitSynch] Added %d synchrotron photons "
            "to photon_list (total now: %d).\n",
            ph_tot, photon_list->num_photons);
    fflush(fPtr);

    return ph_tot;
}
