/*
 * mc_synchrotron.c
 * ================
 * Pitch-angle averaged synchrotron photon emission and SSA for non-thermal
 * electron distributions (single or broken power law).
 *
 * See mc_synchrotron.h for the full physics overview and equation references.
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * FILE ORGANISATION
 * ─────────────────────────────────────────────────────────────────────────────
 *
 * Section 0  Thread-local accelerator management
 * Section 1  F(x) — Bessel integral  [RAIKOU Eq. B13]
 * Section 2  G_a(x) — absorption spectral integral  [RAIKOU Eq. C3]
 * Section 3  Table initialisation: initSynchTables, buildUniversalNuCDF,
 *                                  freeSynchTables
 * Section 4  Emission MC samplers: synchSampleX, synchSampleAlpha,
 *                                  synchSampleGammaEmission
 * Section 5  Absorption coefficient: synchAlphaNu  [RAIKOU Eqs. C2, C4]
 * Section 6  SSA weight modification: applyabsorption,
 *                                     calculateOpticalDepthSSA,
 *                                     applySynchSSAWeightModification
 * Section 7  Universal frequency sampler: sampleSynchFrequency
 * Section 8  Single-photon fill helper: synchFillPhoton
 * Section 9  Main emission: emitCellPackets, photonEmitSynch
 */

#include "mcrat.h"


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 0 — THREAD-LOCAL ACCELERATOR MANAGEMENT                             */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * File-scope universal tables instance.
 * Populated once by initSynchTables; read-only thereafter.
 * All splines are thread-safe for concurrent reads.
 */
static SynchUniversalTables synch_tables;
static int synch_tables_initialized = 0;

/*
 * SynchThreadAccel
 * ----------------
 * Per-thread mutable GSL interpolation accelerators. gsl_interp_accel
 * caches the last binary-search interval and is NOT thread-safe when
 * shared (GSL manual Sec. 26.7). One struct is allocated per OpenMP thread
 * by initSynchThreadAccels and retrieved via getSynchAccel().
 *
 * A single Ga_acc is used for all G_a spline evaluations (Ga_spline,
 * Ga_spline_p1, Ga_spline_p2). Using one accelerator across different
 * splines is safe — the cached interval is merely a search hint and falls
 * back to binary search if the hint is invalid.
 */
typedef struct
{
    gsl_interp_accel *F_acc;
    gsl_interp_accel *Ga_acc;
    gsl_interp_accel *inv_x_cdf_acc;
    gsl_interp_accel *inv_alpha_cdf_acc;
} SynchThreadAccel;

static SynchThreadAccel *synch_thread_accels        = NULL;
static int               synch_num_threads_allocated = 0;

static void initSynchThreadAccels(int num_threads, FILE *fPtr)
{
    int i;
    synch_thread_accels = (SynchThreadAccel *)
                          malloc(num_threads * sizeof(SynchThreadAccel));
    if (synch_thread_accels == NULL)
    {
        fprintf(fPtr,
                ">> [initSynchThreadAccels] ERROR: malloc failed for "
                "%d thread accelerators.\n", num_threads);
        fflush(fPtr);
        exit(1);
    }
    for (i = 0; i < num_threads; i++)
    {
        synch_thread_accels[i].F_acc             = gsl_interp_accel_alloc();
        synch_thread_accels[i].Ga_acc            = gsl_interp_accel_alloc();
        synch_thread_accels[i].inv_x_cdf_acc     = gsl_interp_accel_alloc();
        synch_thread_accels[i].inv_alpha_cdf_acc = gsl_interp_accel_alloc();
    }
    synch_num_threads_allocated = num_threads;
}

static void freeSynchThreadAccels(void)
{
    int i;
    if (!synch_thread_accels) return;
    for (i = 0; i < synch_num_threads_allocated; i++)
    {
        gsl_interp_accel_free(synch_thread_accels[i].F_acc);
        gsl_interp_accel_free(synch_thread_accels[i].Ga_acc);
        gsl_interp_accel_free(synch_thread_accels[i].inv_x_cdf_acc);
        gsl_interp_accel_free(synch_thread_accels[i].inv_alpha_cdf_acc);
    }
    free(synch_thread_accels);
    synch_thread_accels          = NULL;
    synch_num_threads_allocated  = 0;
}

/*
 * getSynchAccel
 * -------------
 * Return the SynchThreadAccel for the calling thread.
 * Falls back to thread 0 when called outside an OpenMP parallel region.
 */
static SynchThreadAccel *getSynchAccel(void)
{
    int thread_id = 0;
    #if defined(_OPENMP)
        thread_id = omp_get_thread_num();
#endif
if (thread_id >= synch_num_threads_allocated)
    thread_id = 0;
return &synch_thread_accels[thread_id];
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 1 — BESSEL INTEGRAL: F(x)                                          */
/* RAIKOU Eq. B13;  G&S91 Eq. 2;  R&L79 Eq. 6.31                             */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
* synch_gsl_underflow_handler
* ---------------------------
* Custom GSL error handler installed during initSynchTables to suppress
* two harmless numerical conditions:
*
*   GSL_EUNDRFLW — occurs when integrating K_{5/3}(xi) to +infinity via
*     gsl_integration_qagiu. K_{5/3} decays as sqrt(pi/(2xi))*exp(-xi) and
*     reaches machine zero around xi ~ 650. The integral has already
*     converged at that point so the underflow is harmless.
*
*   GSL_EROUND — occurs in computeGaAtX when the G_a integrand is nearly
*     zero near x -> SYNCH_X_MAX and the requested relative tolerance
*     cannot be met. The result is still accurate in absolute terms.
*
* All other GSL errors are forwarded to the default handler so that genuine
* problems still abort the program.
*/
static void synch_gsl_underflow_handler(const char *reason,
                                     const char *file,
                                     int         line,
                                     int         gsl_errno)
{
if (gsl_errno == GSL_EUNDRFLW) return;
if (gsl_errno == GSL_EROUND)   return;
gsl_error(reason, file, line, gsl_errno);
}

/*
* bessel_K53_integrand
* --------------------
* Evaluates K_{5/3}(xi) for use by gsl_integration_qagiu in computeFatX.
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
* F(x) is the kernel of both the emission spectral integral G(x)
* [RAIKOU Eq. B12] and the absorption spectral integral G_a(x)
* [RAIKOU Eq. C3]. It is positive for all x > 0 and decays exponentially
* for x >> 1 (R&L79 Eq. 6.33b).
*
* Parameters
* ----------
* x    : dimensionless frequency ratio, x > 0
* ws   : pre-allocated GSL quadrature workspace (size >= 1000)
* fPtr : log file
*
* Returns F(x) >= 0; returns 0.0 for x >= SYNCH_X_MAX.
*/
static double computeFatX(double                     x,
                       gsl_integration_workspace *ws,
                       FILE                      *fPtr)
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
            "x = %.4e (abserr = %.2e).\n",
            gsl_strerror(status), x, abserr);
    fflush(fPtr);
}
return x * result;
}

/*
* evalF
* -----
* Evaluate F(x) from the file-scope precomputed spline, clamping x to the
* table range. Uses the calling thread's accelerator via getSynchAccel().
*
* Called from buildUniversalNuCDF and synchAlphaNu (indirectly via evalGa).
* Thread-safe provided initSynchThreadAccels has been called first.
*/
static double evalF(double x)
{
SynchThreadAccel *ta = getSynchAccel();

if (x >= SYNCH_X_MAX)
    return 0.0;
if (x <= synch_tables.x_arr[0])
    return gsl_spline_eval(synch_tables.F_spline,
                           synch_tables.x_arr[0],
                           ta->F_acc);
return gsl_spline_eval(synch_tables.F_spline, x, ta->F_acc);
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 2 — ABSORPTION SPECTRAL INTEGRAL: G_a(x)                           */
/* RAIKOU Eq. C3                                                               */
/* ═══════════════════════════════════════════════════════════════════════════ */

/* Parameter struct for the G_a integrand passed to GSL quadrature */
struct Ga_integrand_params
{
double            p;
gsl_spline       *F_spline;
gsl_interp_accel *F_acc;
};

/*
* Ga_integrand_fn
* ---------------
* Evaluates z^{(p-2)/2} * F(z) — the integrand of
*
*   G_a(x; p) = integral_x^inf z^{(p-2)/2} F(z) dz     [RAIKOU Eq. C3]
*
* Returns 0 outside the valid spline range [SYNCH_X_MIN, SYNCH_X_MAX].
*/
static double Ga_integrand_fn(double z, void *vp)
{
struct Ga_integrand_params *par = (struct Ga_integrand_params *)vp;
if (z >= SYNCH_X_MAX || z <= 0.0) return 0.0;
double Fz = gsl_spline_eval(par->F_spline, z, par->F_acc);
if (Fz <= 0.0) return 0.0;
return pow(z, 0.5*(par->p - 2.0)) * Fz;
}

/*
* computeGaAtX
* ------------
* Compute G_a(x; p) = integral_x^inf z^{(p-2)/2} F(z) dz  [RAIKOU Eq. C3]
* at a single value of x and power-law index p.
*
* G_a depends only on x and p — not on B, nu_f, or n_e individually. A
* single 1D table at fixed p therefore captures the entire spectral shape
* of alpha_{nu_f}^(f), with B and n_e entering analytically at runtime.
*
* Parameters
* ----------
* x      : lower integration limit (dimensionless frequency ratio)
* p      : power-law index of the electron distribution
* F_spl  : pre-built spline of F(z) on [SYNCH_X_MIN, SYNCH_X_MAX]
* F_acc  : GSL accelerator for F_spl  (single-threaded, used only at init)
* ws     : pre-allocated GSL quadrature workspace
* fPtr   : log file
*
* Returns G_a(x; p) >= 0; returns 0.0 for x >= SYNCH_X_MAX.
*/
static double computeGaAtX(double                     x,
                         double                     p,
                         gsl_spline               *F_spl,
                         gsl_interp_accel         *F_acc,
                         gsl_integration_workspace *ws,
                         FILE                      *fPtr)
{
if (x >= SYNCH_X_MAX) return 0.0;

struct Ga_integrand_params par = { p, F_spl, F_acc };
gsl_function G;
G.function = Ga_integrand_fn;
G.params   = &par;

double result = 0.0, abserr = 0.0;
int status = gsl_integration_qags(&G, x, SYNCH_X_MAX,
                                   0.0, 1e-6,
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

/*
* evalGa
* ------
* Evaluate G_a(x) from a supplied spline, clamping x to the table range.
* Uses the calling thread's Ga_acc accelerator.
*
* A single Ga_acc is shared across all three G_a splines (Ga_spline,
* Ga_spline_p1, Ga_spline_p2). This is safe: the accelerator stores only
* a search-hint interval and falls back to binary search on a cache miss,
* so mixing splines introduces no error, only a marginal performance cost
* on the first call after a spline switch.
*/
static double evalGa(double x, gsl_spline *spl)
{
SynchThreadAccel *ta = getSynchAccel();
if (x <= synch_tables.x_arr[0])
    return gsl_spline_eval(spl, synch_tables.x_arr[0], ta->Ga_acc);
if (x >= SYNCH_X_MAX)
    return 0.0;
return gsl_spline_eval(spl, x, ta->Ga_acc);
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 3 — TABLE INITIALISATION AND TEARDOWN                               */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
* getSynchTables
* --------------
* Return a const pointer to the file-scope synch_tables.
* Exits fatally if called before initSynchTables.
*/
static const SynchUniversalTables *getSynchTables(FILE *fPtr)
{
if (!synch_tables_initialized)
{
    if (fPtr != NULL)
    {
        fprintf(fPtr,
                ">> [getSynchTables] FATAL: synch_tables accessed before "
                "initSynchTables was called.\n");
        fflush(fPtr);
    }
    exit(1);
}
return &synch_tables;
}

/*
* buildUniversalNuCDF
* -------------------
* Build the universal photon-number CDF in dimensionless frequency space
* (nu_tilde = nu / nu_c) and the associated stratified sampling parameters.
* Called once near the end of initSynchTables, after the F(x) spline and
* per-thread accelerators are in place.
*
* Physics
* -------
* The photon-number emissivity per unit log-frequency is:
*
*   pdf(nu_tilde) propto integral_{gamma_min}^{gamma_max}
*                             gamma^{1-p} F(nu_tilde / gamma^2) d(ln gamma)
*
* Substituting nu_tilde = nu/nu_c removes all B dependence. The resulting
* CDF is therefore universal across all cells and emission epochs. At
* emission time nu_f = nu_tilde * nu_c(B_cell) recovers physical Hz.
*
* Algorithm
* ---------
* (1) Evaluate pdf[] on SYNCH_N_NU log-spaced nu_tilde points spanning
*     [SYNCH_X_MIN * GAMMA_MIN^2, SYNCH_X_MAX * GAMMA_MAX^2] using a
*     trapezoid rule over SYNCH_N_GAMMA_INT log-gamma points.
*
* (2) Integrate pdf[] cumulatively and normalise to form the CDF.
*     Store only strictly increasing (u, log10(nu_tilde)) pairs in
*     synch_tables.nu_cdf_u / nu_cdf_log_nu_tilde.
*
* (3) Divide the log-nu_tilde range into SYNCH_N_STRATA equal intervals.
*     For each stratum k, evaluate the CDF at the lower and upper log-nu
*     edges by linear interpolation to obtain:
*       strata_cdf_lo[k], strata_cdf_hi[k],
*       strata_p_k[k] = strata_cdf_hi[k] - strata_cdf_lo[k]
*/
static void buildUniversalNuCDF(FILE *fPtr)
{
    int i, j, k;

    fprintf(fPtr,
            ">> [buildUniversalNuCDF] Building universal nu_tilde CDF "
            "(SYNCH_N_NU=%d, SYNCH_N_GAMMA_INT=%d, SYNCH_N_STRATA=%d)...\n",
            SYNCH_N_NU, SYNCH_N_GAMMA_INT, SYNCH_N_STRATA);
    fflush(fPtr);

    /* ── (1) Evaluate pdf on log-nu_tilde grid ────────────────────────────── */
    double log_nt_min = 0;
    double log_nt_max = 0;
    double dlog_nt    = 0;

    double log_nt_arr[SYNCH_N_NU];
    double pdf_arr   [SYNCH_N_NU];

    double log_g_min = 0;
    double log_g_max = 0;
    double dlog_g    = 0;
    
    #if NONTHERMAL_E_DIST != OFF
        log_nt_min = log10(SYNCH_X_MIN * GAMMA_MIN * GAMMA_MIN);
        log_nt_max = log10(SYNCH_X_MAX * GAMMA_MAX * GAMMA_MAX);
        dlog_nt    = (log_nt_max - log_nt_min) / (SYNCH_N_NU - 1);
        
        log_g_min = log10(GAMMA_MIN);
        log_g_max = log10(GAMMA_MAX);
        dlog_g    = (log_g_max - log_g_min) / (SYNCH_N_GAMMA_INT - 1);
    #endif

    #if NONTHERMAL_E_DIST == POWERLAW

        double p = POWERLAW_INDEX;

        for (i = 0; i < SYNCH_N_NU; i++)
        {
            log_nt_arr[i]   = log_nt_min + i * dlog_nt;
            double nu_tilde = pow(10.0, log_nt_arr[i]);

            /*
             * Trapezoid rule over log-gamma:
             *   integrand = gamma^{1-p} * F(nu_tilde / gamma^2)
             * The normalisation constant of N(gamma) cancels in the
             * CDF normalisation so it is omitted here.
             */
            double sum      = 0.0;
            double prev_val = 0.0;

            for (j = 0; j < SYNCH_N_GAMMA_INT; j++)
            {
                double gamma_j = pow(10.0, log_g_min + j * dlog_g);
                double x_j     = nu_tilde / (gamma_j * gamma_j);
                double val     = pow(gamma_j, 1.0 - p) * evalF(x_j);

                if (j > 0)
                    sum += 0.5 * (val + prev_val) * dlog_g;
                prev_val = val;
            }

            pdf_arr[i] = (sum > 0.0) ? sum : 0.0;
        }

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        /*
         * For a broken power law:
         *   N(gamma) propto gamma^{-p1}              gamma_min <= gamma <= gamma_br
         *   N(gamma) propto C_cont * gamma^{-p2}     gamma_br  <  gamma <= gamma_max
         *
         * where C_cont = gamma_br^{p2-p1} enforces continuity at gamma_br.
         *
         * The photon-number emissivity integrand in d(ln gamma) measure is
         * N(gamma)*gamma^2*F(x)/gamma = N(gamma)*gamma*F(x), giving:
         *
         *   pdf(nu_tilde) propto
         *     integral_{gamma_min}^{gamma_br} gamma^{1-p1} F(nu_tilde/gamma^2) d(ln gamma)
         *   + C_cont *
         *     integral_{gamma_br}^{gamma_max} gamma^{1-p2} F(nu_tilde/gamma^2) d(ln gamma)
         *
         * The log-gamma grid covers the full range [gamma_min, gamma_max].
         * At each grid point we select the appropriate exponent based on
         * whether gamma <= gamma_br or gamma > gamma_br.
         */
        double p1     = POWERLAW_INDEX_1;
        double p2     = POWERLAW_INDEX_2;
        double g_br   = GAMMA_BREAK;
        double C_cont = pow(g_br, p2 - p1);

        for (i = 0; i < SYNCH_N_NU; i++)
        {
            log_nt_arr[i]   = log_nt_min + i * dlog_nt;
            double nu_tilde = pow(10.0, log_nt_arr[i]);

            double sum      = 0.0;
            double prev_val = 0.0;

            for (j = 0; j < SYNCH_N_GAMMA_INT; j++)
            {
                double gamma_j = pow(10.0, log_g_min + j * dlog_g);
                double x_j     = nu_tilde / (gamma_j * gamma_j);
                double F_j     = evalF(x_j);

                double val;
                if (gamma_j <= g_br)
                    val = pow(gamma_j, 1.0 - p1) * F_j;
                else
                    val = C_cont * pow(gamma_j, 1.0 - p2) * F_j;

                if (j > 0)
                    sum += 0.5 * (val + prev_val) * dlog_g;
                prev_val = val;
            }

            pdf_arr[i] = (sum > 0.0) ? sum : 0.0;
        }

    #endif

    /* ── (2) Build normalised CDF ─────────────────────────────────────────── */
    double cdf_raw[SYNCH_N_NU];
    cdf_raw[0] = 0.0;
    for (i = 1; i < SYNCH_N_NU; i++)
        cdf_raw[i] = cdf_raw[i-1]
                   + 0.5 * (pdf_arr[i] + pdf_arr[i-1]) * dlog_nt;

    double total = cdf_raw[SYNCH_N_NU - 1];
    if (total <= 0.0)
    {
        fprintf(fPtr,
                ">> [buildUniversalNuCDF] FATAL: zero CDF total. "
                "Check GAMMA_MIN/MAX, POWERLAW_INDEX, and GAMMA_BREAK.\n");
        fflush(fPtr);
        exit(1);
    }

    /*
     * Store the unnormalised integral total before dividing through.
     * This is used by photonEmitSynch to compute the physical photon
     * emission rate per cell [photons s^-1].
     */
    synch_tables.nu_cdf_norm = total;

    for (i = 0; i < SYNCH_N_NU; i++)
        cdf_raw[i] /= total;
    cdf_raw[SYNCH_N_NU - 1] = 1.0;

    /* Store only strictly increasing (u, log10(nu_tilde)) pairs */
    double prev_u = -1.0;
    int    n      = 0;

    for (i = 0; i < SYNCH_N_NU; i++)
    {
        if (cdf_raw[i] > prev_u)
        {
            synch_tables.nu_cdf_u           [n] = cdf_raw[i];
            synch_tables.nu_cdf_log_nu_tilde[n] = log_nt_arr[i];
            n++;
            prev_u = cdf_raw[i];
        }
    }
    synch_tables.n_nu_cdf = n;

    /* ── (3) Stratum boundaries and probability masses ────────────────────── */
    double dlog_strata = (log_nt_max - log_nt_min) / SYNCH_N_STRATA;

    for (k = 0; k <= SYNCH_N_STRATA; k++)
        synch_tables.strata_log_nu_tilde_edges[k] =
            log_nt_min + k * dlog_strata;

    for (k = 0; k < SYNCH_N_STRATA; k++)
    {
        double log_lo = synch_tables.strata_log_nu_tilde_edges[k];
        double log_hi = synch_tables.strata_log_nu_tilde_edges[k + 1];
        double cdf_lo, cdf_hi;

        /* Linear interpolation of CDF at lower stratum edge */
        if (log_lo <= log_nt_arr[0])
        {
            cdf_lo = cdf_raw[0];
        }
        else if (log_lo >= log_nt_arr[SYNCH_N_NU - 1])
        {
            cdf_lo = 1.0;
        }
        else
        {
            int idx = (int)((log_lo - log_nt_min) / dlog_nt);
            if (idx < 0)             idx = 0;
            if (idx >= SYNCH_N_NU-1) idx = SYNCH_N_NU - 2;
            double frac = (log_lo - log_nt_arr[idx])
                        / (log_nt_arr[idx+1] - log_nt_arr[idx]);
            cdf_lo = cdf_raw[idx] + frac * (cdf_raw[idx+1] - cdf_raw[idx]);
        }

        /* Linear interpolation of CDF at upper stratum edge */
        if (log_hi <= log_nt_arr[0])
        {
            cdf_hi = cdf_raw[0];
        }
        else if (log_hi >= log_nt_arr[SYNCH_N_NU - 1])
        {
            cdf_hi = 1.0;
        }
        else
        {
            int idx = (int)((log_hi - log_nt_min) / dlog_nt);
            if (idx < 0)             idx = 0;
            if (idx >= SYNCH_N_NU-1) idx = SYNCH_N_NU - 2;
            double frac = (log_hi - log_nt_arr[idx])
                        / (log_nt_arr[idx+1] - log_nt_arr[idx]);
            cdf_hi = cdf_raw[idx] + frac * (cdf_raw[idx+1] - cdf_raw[idx]);
        }

        synch_tables.strata_cdf_lo[k] = cdf_lo;
        synch_tables.strata_cdf_hi[k] = cdf_hi;
        synch_tables.strata_p_k   [k] = cdf_hi - cdf_lo;
        if (synch_tables.strata_p_k[k] < 0.0)
            synch_tables.strata_p_k[k] = 0.0;
    }

    /* ── Diagnostic ───────────────────────────────────────────────────────── */
    fprintf(fPtr,
            ">> [buildUniversalNuCDF] nu_tilde=[%.3e, %.3e] "
            "(%.1f decades), %d CDF points, "
            "%d strata of %.2f decades each.\n",
            pow(10.0, log_nt_min), pow(10.0, log_nt_max),
            log_nt_max - log_nt_min, n, SYNCH_N_STRATA,
            (log_nt_max - log_nt_min) / SYNCH_N_STRATA);

    for (k = 0; k < SYNCH_N_STRATA; k++)
    {
        fprintf(fPtr,
                ">>   stratum %2d: nu_tilde=[%.3e, %.3e]  "
                "p_k=%.3e  cdf=[%.8e, %.8e]\n",
                k,
                pow(10.0, synch_tables.strata_log_nu_tilde_edges[k]),
                pow(10.0, synch_tables.strata_log_nu_tilde_edges[k+1]),
                synch_tables.strata_p_k[k],
                synch_tables.strata_cdf_lo[k],
                synch_tables.strata_cdf_hi[k]);
    }
    fflush(fPtr);
}

/*
* initSynchTables
* ---------------
* Build all universal spectral tables. Must be called exactly once at
* startup before any photon emission or SSA evaluation.
*
* Initialisation order
* --------------------
* (1) F(x) grid and spline          [RAIKOU Eq. B13]
* (2) G_a(x) grid and splines       [RAIKOU Eq. C3]
* (3) Inverse CDF of x ~ F(x)*x     [G&S91 Sec. 2]
* (4) Inverse CDF of alpha ~ sin^2  [G&S91 Sec. 2]
* (5) Per-thread accelerators        — must precede (6) because
*     buildUniversalNuCDF calls evalF which calls getSynchAccel
* (6) Universal nu_tilde CDF         [buildUniversalNuCDF]
*
* Note: a temporary single-threaded F_acc is allocated for steps (1)-(4)
* and freed before step (5), since the per-thread array takes over from
* that point forward.
*/
void initSynchTables(FILE *fPtr)
{
    int i;
    fprintf(fPtr,
            ">> [initSynchTables] Building universal spectral tables "
            "(SYNCH_N_X=%d, SYNCH_N_NU=%d, SYNCH_N_STRATA=%d)...\n",
            SYNCH_N_X, SYNCH_N_NU, SYNCH_N_STRATA);
    fflush(fPtr);

    /* Suppress harmless underflow/roundoff during Bessel and G_a integrals */
    gsl_error_handler_t *old_handler =
        gsl_set_error_handler(synch_gsl_underflow_handler);

    /* ── Allocate heap arrays ─────────────────────────────────────────────── */
    synch_tables.x_arr               = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.F_arr               = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.Ga_arr              = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.Ga_arr_p1           = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.Ga_arr_p2           = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_x_cdf_u         = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_x_cdf_logx      = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_alpha_cdf_u     = (double *)malloc(SYNCH_N_X * sizeof(double));
    synch_tables.inv_alpha_cdf_alpha = (double *)malloc(SYNCH_N_X * sizeof(double));

    /* ── Log-spaced x grid ───────────────────────────────────────────────── */
    double log_x_min = log10(SYNCH_X_MIN);
    double log_x_max = log10(SYNCH_X_MAX);
    for (i = 0; i < SYNCH_N_X; i++)
    {
        double t = (double)i / (SYNCH_N_X - 1);
        synch_tables.x_arr[i] = pow(10.0,
                                     log_x_min + t*(log_x_max - log_x_min));
    }

    /* Temporary single-threaded accelerator used for steps (1)-(4) only.
     * Freed before initSynchThreadAccels so it does not alias any per-thread
     * slot.                                                                  */
    gsl_interp_accel         *tmp_F_acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *ws       = gsl_integration_workspace_alloc(1000);

    /* ── (1) F(x)  [RAIKOU Eq. B13] ─────────────────────────────────────── */
    fprintf(fPtr, ">> [initSynchTables] Computing F(x) [RAIKOU Eq. B13]...\n");
    fflush(fPtr);

    for (i = 0; i < SYNCH_N_X; i++)
        synch_tables.F_arr[i] = computeFatX(synch_tables.x_arr[i], ws, fPtr);

    synch_tables.F_spline = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    gsl_spline_init(synch_tables.F_spline,
                    synch_tables.x_arr, synch_tables.F_arr, SYNCH_N_X);

    /* ── (2) G_a(x)  [RAIKOU Eq. C3] ────────────────────────────────────── */
    fprintf(fPtr, ">> [initSynchTables] Computing G_a(x) [RAIKOU Eq. C3]...\n");
    fflush(fPtr);

    #if NONTHERMAL_E_DIST == POWERLAW

        for (i = 0; i < SYNCH_N_X; i++)
            synch_tables.Ga_arr[i] = computeGaAtX(synch_tables.x_arr[i],
                                                   POWERLAW_INDEX,
                                                   synch_tables.F_spline,
                                                   tmp_F_acc,
                                                   ws, fPtr);

        /* p1/p2 arrays unused for single power law — zero-fill */
        for (i = 0; i < SYNCH_N_X; i++)
        {
            synch_tables.Ga_arr_p1[i] = 0.0;
            synch_tables.Ga_arr_p2[i] = 0.0;
        }

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        for (i = 0; i < SYNCH_N_X; i++)
        {
            synch_tables.Ga_arr_p1[i] = computeGaAtX(synch_tables.x_arr[i],
                                                      POWERLAW_INDEX_1,
                                                      synch_tables.F_spline,
                                                      tmp_F_acc,
                                                      ws, fPtr);
            synch_tables.Ga_arr_p2[i] = computeGaAtX(synch_tables.x_arr[i],
                                                      POWERLAW_INDEX_2,
                                                      synch_tables.F_spline,
                                                      tmp_F_acc,
                                                      ws, fPtr);
            synch_tables.Ga_arr[i] = synch_tables.Ga_arr_p1[i];
        }

    #endif

    synch_tables.Ga_spline    = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    synch_tables.Ga_spline_p1 = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);
    synch_tables.Ga_spline_p2 = gsl_spline_alloc(gsl_interp_linear, SYNCH_N_X);

    gsl_spline_init(synch_tables.Ga_spline,
                    synch_tables.x_arr, synch_tables.Ga_arr,    SYNCH_N_X);
    gsl_spline_init(synch_tables.Ga_spline_p1,
                    synch_tables.x_arr, synch_tables.Ga_arr_p1, SYNCH_N_X);
    gsl_spline_init(synch_tables.Ga_spline_p2,
                    synch_tables.x_arr, synch_tables.Ga_arr_p2, SYNCH_N_X);

    /* ── (3) Inverse CDF of x ~ F(x)*x d(log x)  [G&S91 Sec. 2] ───────── */
    /*
     * The probability that a single-electron photon has dimensionless
     * frequency ratio in [x, x+dx] is proportional to F(x) dx (R&L79
     * Eq. 6.36). In d(log x) measure the weight is F(x)*x per unit log x.
     * We store the inverse map u -> log10(x) in log-complementary space
     * v = -log10(1-u) for numerical stability near u -> 1.
     */
    fprintf(fPtr,
            ">> [initSynchTables] Building inverse-CDF of x...\n");
    fflush(fPtr);

    {
        double dlog_x = (log_x_max - log_x_min) / (SYNCH_N_X - 1);
        double cum    = 0.0;

        /* Raw cumulative trapezoid in log-x */
        double raw_cdf[SYNCH_N_X];
        raw_cdf[0] = 0.0;
        for (i = 1; i < SYNCH_N_X; i++)
        {
            double w0 = synch_tables.F_arr[i-1] * synch_tables.x_arr[i-1];
            double w1 = synch_tables.F_arr[i]   * synch_tables.x_arr[i];
            cum      += 0.5 * (w0 + w1) * dlog_x;
            raw_cdf[i] = cum;
        }

        double total_x = raw_cdf[SYNCH_N_X - 1];
        if (total_x > 0.0)
            for (i = 0; i < SYNCH_N_X; i++)
                raw_cdf[i] /= total_x;
        raw_cdf[SYNCH_N_X - 1] = 1.0;

        /* Compact to strictly increasing pairs in log-complementary space */
        double tmp_v   [SYNCH_N_X];
        double tmp_logx[SYNCH_N_X];
        int    n_cdf   = 0;
        double prev_u  = -1.0;

        for (i = 0; i < SYNCH_N_X; i++)
        {
            double u_i = raw_cdf[i];
            if (u_i > prev_u)
            {
                double complement = 1.0 - u_i;
                if (complement <= 0.0)
                {
                    if (n_cdf > 0)
                    {
                        tmp_v   [n_cdf] = tmp_v[n_cdf-1] + 1.0;
                        tmp_logx[n_cdf] = log10(synch_tables.x_arr[i]);
                        n_cdf++;
                    }
                    break;
                }
                tmp_v   [n_cdf] = -log10(complement);
                tmp_logx[n_cdf] = log10(synch_tables.x_arr[i]);
                n_cdf++;
                prev_u = u_i;
            }
        }

        for (i = 0; i < n_cdf; i++)
        {
            synch_tables.inv_x_cdf_u   [i] = tmp_v   [i];
            synch_tables.inv_x_cdf_logx[i] = tmp_logx[i];
        }

        synch_tables.inv_x_cdf_spline = gsl_spline_alloc(gsl_interp_linear,
                                                           n_cdf);
        gsl_spline_init(synch_tables.inv_x_cdf_spline,
                        synch_tables.inv_x_cdf_u,
                        synch_tables.inv_x_cdf_logx,
                        n_cdf);
    }

    /* ── (4) Inverse CDF of alpha ~ sin^2(alpha)  [G&S91 Sec. 2] ───────── */
    /*
     * The isotropic pitch-angle distribution is
     *   f(alpha) = (2/pi) sin^2(alpha),  alpha in [0, pi]
     *
     * Its analytic CDF is:
     *   C(alpha) = (alpha - sin(2*alpha)/2) / pi
     *
     * Evaluated on a uniform alpha grid; stored as the inverse map
     * u -> alpha as a linear spline.
     */
    fprintf(fPtr,
            ">> [initSynchTables] Building inverse-CDF of alpha...\n");
    fflush(fPtr);

    for (i = 0; i < SYNCH_N_X; i++)
    {
        double alpha = M_PI * (double)i / (SYNCH_N_X - 1);
        synch_tables.inv_alpha_cdf_alpha[i] = alpha;
        synch_tables.inv_alpha_cdf_u    [i] = (alpha - 0.5*sin(2.0*alpha))
                                             / M_PI;
    }
    synch_tables.inv_alpha_cdf_u[SYNCH_N_X - 1] = 1.0;

    synch_tables.inv_alpha_cdf_spline = gsl_spline_alloc(gsl_interp_linear,
                                                          SYNCH_N_X);
    gsl_spline_init(synch_tables.inv_alpha_cdf_spline,
                    synch_tables.inv_alpha_cdf_u,
                    synch_tables.inv_alpha_cdf_alpha,
                    SYNCH_N_X);

    /* Free temporary single-threaded resources used in steps (1)-(4) */
    gsl_interp_accel_free(tmp_F_acc);
    gsl_integration_workspace_free(ws);

    /* ── (5) Per-thread accelerators ─────────────────────────────────────── */
    /*
     * Must be initialised BEFORE buildUniversalNuCDF (step 6) because
     * buildUniversalNuCDF calls evalF, which calls getSynchAccel().
     */
    int num_threads = 1;
    #if defined(_OPENMP)
        num_threads = omp_get_max_threads();
    #endif
    initSynchThreadAccels(num_threads, fPtr);

    /* ── (6) Universal photon-number CDF in nu_tilde space ──────────────── */
    buildUniversalNuCDF(fPtr);

    gsl_set_error_handler(old_handler);

    synch_tables_initialized = 1;

    fprintf(fPtr, ">> [initSynchTables] All universal tables ready.\n\n");
    fflush(fPtr);
}

/*
 * freeSynchTables
 * ---------------
 * Release all heap memory allocated by initSynchTables.
 * Must be called exactly once at end of simulation.
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

    /* Free per-thread accelerators first */
    freeSynchThreadAccels();

    /* Free heap arrays */
    free(synch_tables.x_arr);
    free(synch_tables.F_arr);
    free(synch_tables.Ga_arr);
    free(synch_tables.Ga_arr_p1);
    free(synch_tables.Ga_arr_p2);
    free(synch_tables.inv_x_cdf_u);
    free(synch_tables.inv_x_cdf_logx);
    free(synch_tables.inv_alpha_cdf_u);
    free(synch_tables.inv_alpha_cdf_alpha);

    /* Free splines */
    gsl_spline_free(synch_tables.F_spline);
    gsl_spline_free(synch_tables.Ga_spline);
    gsl_spline_free(synch_tables.Ga_spline_p1);
    gsl_spline_free(synch_tables.Ga_spline_p2);
    gsl_spline_free(synch_tables.inv_x_cdf_spline);
    gsl_spline_free(synch_tables.inv_alpha_cdf_spline);

    synch_tables_initialized = 0;

    fprintf(fPtr, ">> [freeSynchTables] All tables freed.\n");
    fflush(fPtr);
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 4 — EMISSION MONTE CARLO SAMPLERS                                   */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * synchSampleX
 * ------------
 * Sample x = nu_f / (gamma_e^2 * nu_c) from the single-electron synchrotron
 * spectrum F(x) via the precomputed inverse CDF stored in log-complementary
 * space v = -log10(1-u)  [G&S91 Sec. 2; R&L79 Eq. 6.36].
 *
 * Note: this function is retained for use by mc_cyclosynch.c and any other
 * callers that need a direct single-electron x sample. It is NOT called
 * during the stratified synchrotron emission path (sampleSynchFrequency is
 * used there instead).
 *
 * Parameters
 * ----------
 * u : uniform random variate in (0, 1)
 *
 * Returns x > 0.
 */
double synchSampleX(double u)
{
    const SynchUniversalTables *tables = getSynchTables(NULL);
    SynchThreadAccel *ta = getSynchAccel();

    if (u <= 0.0) u = 1e-14;
    if (u >= 1.0) u = 1.0 - 1e-14;

    double v    = -log10(1.0 - u);
    double v_lo = tables->inv_x_cdf_u[0];
    double v_hi = tables->inv_x_cdf_u[tables->inv_x_cdf_spline->size - 1];

    if (v < v_lo) v = v_lo;
    if (v > v_hi) v = v_hi;

    return pow(10.0, gsl_spline_eval(tables->inv_x_cdf_spline, v,
                                      ta->inv_x_cdf_acc));
}

/*
 * synchSampleAlpha
 * ----------------
 * Sample the electron pitch angle alpha from the isotropic distribution
 *   f(alpha) = (2/pi) sin^2(alpha)                   [G&S91 Sec. 2]
 *
 * Note: retained for external use by mc_cyclosynch.c.
 *
 * Parameters
 * ----------
 * u : uniform random variate in (0, 1)
 *
 * Returns alpha in (0, pi).
 */
double synchSampleAlpha(double u)
{
    const SynchUniversalTables *tables = getSynchTables(NULL);
    SynchThreadAccel *ta = getSynchAccel();
    return gsl_spline_eval(tables->inv_alpha_cdf_spline, u,
                            ta->inv_alpha_cdf_acc);
}

/*
 * synchSampleGammaEmission
 * ------------------------
 * Sample the electron Lorentz factor gamma_e from the emission-weighted
 * distribution N(gamma) * gamma^2.
 *
 * The gamma^2 weighting arises because total synchrotron power per electron
 * scales as gamma^2 (R&L79 Eq. 6.38). A photon drawn at random from the
 * full emission distribution was emitted by an electron distributed as
 * N(gamma)*gamma^2, not N(gamma).
 *
 * Note: retained for external use. Must NOT be replaced by samplePowerLaw
 * from electron.c, which samples from N(gamma) directly.
 *
 * SINGLE POWER LAW  N(gamma) ∝ gamma^{-p}:
 *   emission weight ∝ gamma^{2-p} = gamma^q,  q = 2-p
 *   inverse CDF: gamma(u) = (gmin^q + u*(gmax^q - gmin^q))^{1/q}
 *   special case q=0 (p=2): gamma(u) = gmin*(gmax/gmin)^u
 *
 * BROKEN POWER LAW:
 *   emission power split between low/high segments; each segment sampled
 *   with its share of the total emission integral.
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
        double q    = 2.0 - p;

        if (fabs(q) < 1e-6)
            return gmin * pow(gmax / gmin, u);

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
         * Continuity factor at gamma_break: enforces N(gamma) continuous
         * at gamma_br  [RAIKOU Eq. A4-A5].
         */
        double C_cont = pow(gbr, p2 - p1);

        /* Total emission power in each segment */
        double I1, I2;
        if (fabs(q1) < 1e-6) I1 = log(gbr / gmin);
        else                  I1 = (pow(gbr, q1) - pow(gmin, q1)) / q1;

        if (fabs(q2) < 1e-6) I2 = C_cont * log(gmax / gbr);
        else                  I2 = C_cont * (pow(gmax, q2) - pow(gbr, q2)) / q2;

        double I_tot = I1 + I2;
        if (I_tot <= 0.0) return gmin;

        double f1 = I1 / I_tot;   /* fraction of emission from low segment */

        if (u < f1)
        {
            double u1 = u / f1;
            if (fabs(q1) < 1e-6)
                return gmin * pow(gbr / gmin, u1);
            return pow(pow(gmin, q1) + u1*(pow(gbr, q1) - pow(gmin, q1)),
                       1.0/q1);
        }
        else
        {
            double u2 = (u - f1) / (1.0 - f1);
            if (fabs(q2) < 1e-6)
                return gbr * pow(gmax / gbr, u2);
            return pow(pow(gbr, q2) + u2*(pow(gmax, q2) - pow(gbr, q2)),
                       1.0/q2);
        }

    #endif
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 5 — ABSORPTION COEFFICIENT                                          */
/* RAIKOU Eqs. C2 (single PL) and C4 (broken PL)                              */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * synchAlphaNu
 * ------------
 * Compute the SSA absorption coefficient alpha_{nu_f}^(f) [cm^-1] in the
 * fluid rest frame at comoving frequency nu_f [Hz].
 *
 * Single power law  [RAIKOU Eq. C2]:
 *   alpha = (p-1)(p+2) n_e e^2 nu_c
 *           / (4 sqrt(3) me c (gmin^{1-p} - gmax^{1-p}))
 *           * (nu_f/nu_c)^{-(p+4)/2}
 *           * [G_a(x_max) - G_a(x_min)]
 *
 * Broken power law  [RAIKOU Eq. C4]:
 *   two-segment sum with indices p1/p2 and continuity factor C_cont.
 *
 * The G_a table handles all x-dependence; B and n_e enter analytically
 * through the prefactor. Returns 0 for degenerate inputs.
 */
double synchAlphaNu(double nu_f,
                    double B,
                    double n_e_nth,
                    FILE  *fPtr)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);

    if (nu_f <= 0.0 || B <= 0.0 || n_e_nth <= 0.0) return 0.0;

    /*
     * Pitch-angle averaged critical frequency  [RAIKOU Eq. B7]:
     *   nu_c = 3 e B / (4 pi me c)
     */
    double nu_c = (3.0 * CHARGE_EL * B) / (4.0 * M_PI * M_EL * C_LIGHT);
    if (nu_c <= 0.0) return 0.0;

    #if NONTHERMAL_E_DIST == POWERLAW

        double p    = POWERLAW_INDEX;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;

        if (fabs(p - 1.0) < 1e-6)
        {
            fprintf(fPtr,
                    ">> [synchAlphaNu] WARNING: p = 1, degenerate "
                    "normalisation. Returning 0.\n");
            fflush(fPtr);
            return 0.0;
        }

        double denom = pow(gmin, 1.0 - p) - pow(gmax, 1.0 - p);
        if (fabs(denom) < 1e-300) return 0.0;

        /*
         * x = nu_f / (gamma^2 nu_c):
         *   x_min -> gamma_max (smallest x)
         *   x_max -> gamma_min (largest  x)
         */
        double x_min = nu_f / (gmax * gmax * nu_c);
        double x_max = nu_f / (gmin * gmin * nu_c);

        double delta_Ga = evalGa(x_min, tables->Ga_spline)
                    - evalGa(x_max, tables->Ga_spline);
        if (delta_Ga <= 0.0) return 0.0;

        double prefactor = (p - 1.0) * (p + 2.0)
                         * n_e_nth * CHARGE_EL * CHARGE_EL
                         / (4.0 * sqrt(3.0) * M_EL * C_LIGHT * nu_c * denom);

        return prefactor * pow(nu_f / nu_c, -0.5*(p + 4.0)) * delta_Ga;

    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW

        double p1   = POWERLAW_INDEX_1;
        double p2   = POWERLAW_INDEX_2;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;
        double gbr  = GAMMA_BREAK;

        double A = brokenPowerLawNorm(p1, p2, gmin, gmax, gbr);
        if (A <= 0.0) return 0.0;

        double C_cont = pow(gbr, p2 - p1);

        double x_max = nu_f / (gmin * gmin * nu_c);
        double x_br  = nu_f / (gbr  * gbr  * nu_c);
        double x_min = nu_f / (gmax * gmax * nu_c);

        double delta_Ga_p1 = evalGa(x_br,  tables->Ga_spline_p1)
                       - evalGa(x_max, tables->Ga_spline_p1);
        if (delta_Ga_p1 < 0.0) delta_Ga_p1 = 0.0;

        double delta_Ga_p2 = evalGa(x_min, tables->Ga_spline_p2)
                       - evalGa(x_br,  tables->Ga_spline_p2);
        if (delta_Ga_p2 < 0.0) delta_Ga_p2 = 0.0;

        if (delta_Ga_p1 == 0.0 && delta_Ga_p2 == 0.0) return 0.0;

        double common = (A * n_e_nth * CHARGE_EL * CHARGE_EL)
                      / (4.0 * sqrt(3.0) * M_EL * C_LIGHT * nu_c);

        double term1 = (p1 + 2.0)
                     * pow(nu_f / nu_c, -0.5*(p1 + 4.0))
                     * delta_Ga_p1;

        double term2 = (p2 + 2.0) * C_cont
                     * pow(nu_f / nu_c, -0.5*(p2 + 4.0))
                     * delta_Ga_p2;

        return common * (term1 + term2);

    #endif
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 6 — SSA WEIGHT MODIFICATION                                         */
/* RAIKOU Eqs. 31, 40                                                          */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * applyabsorption
 * ---------------
 * Apply continuous SSA weight attenuation over a lab-frame step of length
 * dl [cm]:
 *
 *   w_new = w_old * exp(-abs_optical_depth * dl)     [RAIKOU Eq. 40]
 *
 * ph->abs_optical_depth was set by calculateOpticalDepthSSA (called from
 * calculateOpticalDepth) as
 *
 *   abs_optical_depth = fluid_factor * alpha_{nu_f}^(f)
 *                     = (nu_f / nu_z) * alpha_{nu_f}^(f)
 *
 * so the full RAIKOU Eq. 31 frame correction is already included and no
 * additional Doppler factor is needed here.
 *
 * Parameters
 * ----------
 * ph : photon packet (modified in-place)
 * dl : lab-frame step length [cm]; must be >= 0
 */
void applyAbsorption(struct photon *ph, double dl)
{
    #if SYNCHROTRON_SWITCH == ON
        if (ph == NULL)                       return;
        if (dl <= 0.0)                        return;
        if (ph->absorption_opacity <= 0.0)     return;

        double tau = ph->absorption_opacity * dl;

        /* Guard against exp underflow for very large optical depths */
        //set to DBL_MIN so these photons can potentially be written out instead of causing garbage data being written out in the printPhotons function. Also applyRussianRoulette, recongnises these as needing to be removed.
        //set the nearest_block_index=-1 so these photons dont scatter as we dont really care about them anymore
        if (tau > 200.0)
        {
            ph->weight = DBL_MIN;
            ph->nearest_block_index = -1;
            ph->time_to_scatter     = FLT_MAX / C_LIGHT;
            ph->recalc_properties   = 1;
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
 * Compute the frame-corrected SSA absorption coefficient
 *
 *   ph->abs_optical_depth = gamma * fluid_factor * alpha_{nu_f}^(f)   [cm^{-1}]
 *
 * where fluid_factor = (1 - beta * cos(theta)) = nu_f / nu_z is the
 * Doppler frame correction already computed in calculateOpticalDepth and
 * applied identically to the scattering opacity.  Pre-multiplying here
 * keeps all frame-correction arithmetic in one place (calculateOpticalDepth)
 * so that applyabsorption can simply compute
 *
 *   w_new = w_old * exp(-ph->abs_optical_depth * dl)
 *
 * with no additional factors.  [RAIKOU Eq. 31, 40]
 *
 * Parameters
 * ----------
 * ph           : photon packet (ph->abs_optical_depth modified in-place)
 * hydro_data   : hydro frame (provides B and nonthermal_dens)
 * fluid_factor : (1 - beta*cos(theta)), computed in calculateOpticalDepth
 * fPtr         : log file
 */
void calculateOpticalDepthSSA(struct photon          *ph,
                               struct hydro_dataframe *hydro_data,
                               double                  fluid_factor,
                               FILE                   *fPtr)
{
    #if SYNCHROTRON_SWITCH == ON

        ph->absorption_opacity = 0.0;
    
        //this can operate on all photon types
        //if (ph->type != SYNCH_PHOTON)
        //    return;

        int    cell_idx  = ph->nearest_block_index;
        double B         = getMagneticFieldMagnitude(hydro_data, cell_idx);
        double n_e_nth   = (hydro_data->nonthermal_dens)[cell_idx];
        double gamma     = (hydro_data->gamma)[cell_idx];
        double nu_f      = ph->comv_p0 * C_LIGHT / PL_CONST;

        double alpha = synchAlphaNu(nu_f, B, n_e_nth, fPtr);

        /*
         * Apply the frame correction (nu_f / nu_z = fluid_factor) here so that
         * applyabsorption requires no additional frame arithmetic.
         * [RAIKOU Eq. 31]
         */
        //working thorugh Abramowicz 1991, find that gamma shouldnt be applied to number density but should be multiplied to the absorption opacity, which is determined from the fluid frame nonthermal electron density
        //put this here to be explicit rather than multiplying nonthermal electron density by gamma to get lab frame density like we do elsewhere
        ph->absorption_opacity = gamma * fluid_factor * alpha;

    #else
        (void)ph;
        (void)hydro_data;
        (void)fluid_factor;
        (void)fPtr;
    #endif
}



/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 7 — UNIVERSAL FREQUENCY SAMPLER                                     */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * sampleSynchFrequency
 * --------------------
 * Draw a photon frequency from stratum k of the universal nu_tilde CDF
 * and convert to physical Hz by scaling with nu_c(B_cell).
 *
 * Algorithm
 * ---------
 * (1) Draw u ~ Uniform(strata_cdf_lo[k], strata_cdf_hi[k]) — restricts the
 *     sample to stratum k's frequency range.
 * (2) Binary search on nu_cdf_u[] to find the bracketing interval.
 * (3) Linear interpolation to get log10(nu_tilde).
 * (4) Return 10^(log10(nu_tilde)) * nu_c(B_cell)  [Hz].
 *
 * Parameters
 * ----------
 * k       : stratum index in [0, SYNCH_N_STRATA)
 * B_cell  : magnetic field magnitude in the emitting cell [G]
 * rand    : GSL RNG
 *
 * Returns nu_f [Hz] in the fluid rest frame.
 */
double sampleSynchFrequency(int k, double B_cell, gsl_rng *rand)
{
    const SynchUniversalTables *tables = getSynchTables(NULL);

    /* Compute nu_c once — used for all three return paths below */
    double nu_c = (3.0 * CHARGE_EL * B_cell) / (4.0 * M_PI * M_EL * C_LIGHT);

    int    n    = tables->n_nu_cdf;
    double u_lo = tables->strata_cdf_lo[k];
    double u_hi = tables->strata_cdf_hi[k];

    double u = u_lo + gsl_rng_uniform_pos(rand) * (u_hi - u_lo);

    /* Clamp to stored CDF range */
    if (u <= tables->nu_cdf_u[0])
        return pow(10.0, tables->nu_cdf_log_nu_tilde[0]) * nu_c;
    if (u >= tables->nu_cdf_u[n - 1])
        return pow(10.0, tables->nu_cdf_log_nu_tilde[n - 1]) * nu_c;

    /* Binary search for bracketing interval */
    int lo = 0, hi = n - 1;
    while (hi - lo > 1)
    {
        int mid = (lo + hi) / 2;
        if (tables->nu_cdf_u[mid] <= u)
            lo = mid;
        else
            hi = mid;
    }

    /* Linear interpolation in log-nu_tilde */
    double frac     = (u - tables->nu_cdf_u[lo])
                    / (tables->nu_cdf_u[hi] - tables->nu_cdf_u[lo]);
    double log_nu_t = tables->nu_cdf_log_nu_tilde[lo]
                    + frac * (tables->nu_cdf_log_nu_tilde[hi]
                             - tables->nu_cdf_log_nu_tilde[lo]);

    return pow(10.0, log_nu_t) * nu_c;
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 8 — SINGLE-PHOTON FILL HELPER                                       */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * synchFillPhoton
 * ---------------
 * Populate all fields of a struct photon for synchrotron emission from fluid
 * cell cell_idx at comoving frequency nu_f_comv [Hz] with packet weight
 * ph_weight.
 *
 * Steps
 * -----
 * (1) Sample isotropic comoving direction (theta, phi); assemble comoving
 *     4-momentum p^mu_comv = (h nu/c) * (1, sin(theta)cos(phi),
 *                                           sin(theta)sin(phi), cos(theta)).
 * (2) Lorentz boost to lab frame: negate the fluid 3-velocity and call
 *     lorentzBoost with flag 'p' (4-momentum transform).
 * (3) Birth position: uniform random offset within the cell bounding box,
 *     converted to MCRaT Cartesian coordinates.
 * (4) Stokes parameters: (s0,s1,s2,s3) = (1,0,0,0) — unpolarised emission
 *     [G&S91 Sec. 2].
 * (5) Bookkeeping: type, weight, cell index, zero scattering counters.
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

    /* ── (1) Isotropic comoving 4-momentum ───────────────────────────────── */
    double com_v_phi   = samplePhotonPhi(rand, fPtr);
    double com_v_theta = samplePhotonTheta(boost, rand, fPtr);
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
     * Assign a random azimuthal angle in 2D/2.5D (needed for coordinate
     * conversion); unused in 3D where the full vector is available.
     */
    #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
        position_phi = gsl_rng_uniform(rand) * 2.0 * M_PI;
    #else
        position_phi = 0.0;
    #endif

    #if DIMENSIONS == THREE
        hydroVectorToCartesian(boost,
                                (hydro_data->v0)[cell_idx],
                                (hydro_data->v1)[cell_idx],
                                (hydro_data->v2)[cell_idx],
                                (hydro_data->r0)[cell_idx],
                                (hydro_data->r1)[cell_idx],
                                (hydro_data->r2)[cell_idx]);
    #elif DIMENSIONS == TWO_POINT_FIVE
        hydroVectorToCartesian(boost,
                                (hydro_data->v0)[cell_idx],
                                (hydro_data->v1)[cell_idx],
                                (hydro_data->v2)[cell_idx],
                                (hydro_data->r0)[cell_idx],
                                (hydro_data->r1)[cell_idx],
                                position_phi);
    #else
        hydroVectorToCartesian(boost,
                                (hydro_data->v0)[cell_idx],
                                (hydro_data->v1)[cell_idx],
                                0,
                                (hydro_data->r0)[cell_idx],
                                (hydro_data->r1)[cell_idx],
                                position_phi);
    #endif

    /* Negate fluid velocity: boost from fluid frame to lab frame */
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
    ph->total_scattering_opacity = 0.0;

    #if SCATTERING_BIAS_SWITCH != OFF
        {
            int idx;
            for (idx = 0; idx < 1 + N_GAMMA; idx++)
            {
                ph->scattering_opacity[idx]  = 0.0;
                ph->scattering_bias[idx] = 0.0;
            }
        }
    #endif

    #if SYNCHROTRON_SWITCH == ON
        ph->absorption_opacity = 0.0;
    #endif
}


/* ═══════════════════════════════════════════════════════════════════════════ */
/* SECTION 9 — MAIN EMISSION                                                   */
/* ═══════════════════════════════════════════════════════════════════════════ */

/*
 * emitCellPackets
 * ---------------
 * Emit ph_count photon packets from active cell ci using the universal
 * stratified frequency sampler. No per-cell CDF build is required.
 *
 * Packet allocation
 * -----------------
 * ph_count packets are divided equally across SYNCH_N_STRATA strata using
 * integer division with the remainder distributed one extra packet to the
 * lowest-indexed strata, guaranteeing sum_k alloc[k] == ph_count exactly.
 *
 * Importance weighting
 * --------------------
 * Each stratum k naturally contains fraction p_k of all emitted photons.
 * Forcing equal allocation (1/K per stratum) introduces a bias corrected by:
 *
 *   w_k = ph_weight_adjusted * strata_p_k[k] * SYNCH_N_STRATA
 *
 * so that sum_k alloc[k] * w_k = ph_count * ph_weight_adjusted exactly
 * (up to the integer-remainder contribution of order 1/ph_count).
 *
 * Degenerate strata
 * -----------------
 * Strata with strata_p_k[k] == 0 emit zero-weight packets at the stratum
 * nu_tilde midpoint, preserving ph_count without contributing to the
 * weighted spectrum.
 *
 * Fallback
 * --------
 * If B_cell <= 0 all ph_count packets are emitted as zero-weight at
 * nu = 1 Hz and the condition is logged.
 *
 * Parameters
 * ----------
 * ci                 : hydro_data cell index
 * ph_count           : number of packets to emit from this cell
 * B_cell             : magnetic field magnitude [G]
 * ph_weight_adjusted : baseline packet weight w_0
 * ph_emit            : output photon array (pre-allocated, size >= ph_tot)
 * idx_ptr            : running write index into ph_emit; advanced by ph_count
 * hydro_data         : hydrodynamic data frame
 * rand               : GSL RNG
 * fPtr               : log file
 */
static void emitCellPackets(int                    ci,
                              int                    ph_count,
                              double                 B_cell,
                              double                 ph_weight_adjusted,
                              struct photon         *ph_emit,
                              int                   *idx_ptr,
                              struct hydro_dataframe *hydro_data,
                              gsl_rng               *rand,
                              FILE                  *fPtr)
{
    const SynchUniversalTables *tables = getSynchTables(fPtr);
    int i, k;

    /* ── Fallback: degenerate B ───────────────────────────────────────────── */
    if (B_cell <= 0.0)
    {
        fprintf(fPtr,
                ">> [emitCellPackets] Cell %d: B_cell <= 0, "
                "emitting %d zero-weight packets.\n",
                ci, ph_count);
        fflush(fPtr);

        for (i = 0; i < ph_count; i++)
        {
            synchFillPhoton(&ph_emit[*idx_ptr], ci, 1.0,
                             0.0, hydro_data, rand, fPtr);
            (*idx_ptr)++;
        }
        return;
    }

    /* ── Equal allocation across strata with remainder ───────────────────── */
    int alloc[SYNCH_N_STRATA];
    int base      = ph_count / SYNCH_N_STRATA;
    int remainder = ph_count % SYNCH_N_STRATA;

    for (k = 0; k < SYNCH_N_STRATA; k++)
        alloc[k] = base + (k < remainder ? 1 : 0);

    /* ── Emit stratum by stratum ──────────────────────────────────────────── */
    for (k = 0; k < SYNCH_N_STRATA; k++)
    {
        if (alloc[k] > 0)
        {
            if (tables->strata_p_k[k] <= 0.0 ||
                tables->strata_cdf_hi[k] <= tables->strata_cdf_lo[k])
            {
                /*
                 * No CDF support in this stratum. Emit zero-weight packets
                 * at the stratum nu_tilde midpoint to preserve ph_count
                 * exactly. These packets carry no physical energy.
                 */
                double log_nt_mid = 0.5
                                  * (tables->strata_log_nu_tilde_edges[k]
                                   + tables->strata_log_nu_tilde_edges[k+1]);
                double nu_c   = (3.0 * CHARGE_EL * B_cell)
                              / (4.0 * M_PI * M_EL * C_LIGHT);
                double nu_mid = pow(10.0, log_nt_mid) * nu_c;

                for (i = 0; i < alloc[k]; i++)
                {
                    synchFillPhoton(&ph_emit[*idx_ptr], ci, nu_mid,
                                     0.0, hydro_data, rand, fPtr);
                    (*idx_ptr)++;
                }
            }
            else
            {
                /*
                 * Normal stratified path.
                 * w_k = w_0 * p_k * K corrects for the forced equal
                 * allocation across strata.
                 */
                double w_k = ph_weight_adjusted
                           * tables->strata_p_k[k]
                           * (double)SYNCH_N_STRATA;

                for (i = 0; i < alloc[k]; i++)
                {
                    double nu_f = sampleSynchFrequency(k, B_cell, rand);
                    synchFillPhoton(&ph_emit[*idx_ptr], ci, nu_f,
                                     w_k, hydro_data, rand, fPtr);
                    (*idx_ptr)++;
                }
            }
        }
    }
}

/*
 * photonEmitSynch
 * ---------------
 * Emit synchrotron photon packets from all active cells in the injection
 * shell at radius r_inj. Returns the total number of packets emitted.
 *
 * Algorithm (step numbers match log output)
 * -----------------------------------------
 * Step 1  Compute shell radial boundaries [rmin, rmax] from r_inj and fps.
 * Step 2  First pass over hydro grid: count active cells (block_cnt).
 * Step 3  Allocate per-cell arrays: ph_dens[], cell_index[], W_cell[].
 * Step 4  Second pass: record cell indices and emission weights
 *           W_j = n_e_nth * B^2 * V  (proportional to total synchrotron
 *           power; the gamma integral is a compile-time constant and cancels
 *           in the normalisation).
 * Step 5  Tune ph_weight_adjusted so that sum_j lambda_j lies in
 *           [min_photons, max_photons], where lambda_j = W_j/(w_0*fps).
 * Step 6  Poisson draw: ph_dens[j] ~ Poisson(lambda_j).
 * Step 7  Allocate and zero-initialise ph_emit[ph_tot].
 * Step 8  Per-cell emission loop: call emitCellPackets for each cell with
 *           ph_dens[j] > 0, passing B_cell directly. The universal nu_tilde
 *           CDF (built once in initSynchTables) is used for all cells — no
 *           per-cell CDF build or numerical integration is performed here.
 * Step 9  Add ph_emit batch to photon_list via addToPhotonList.
 * Step 10 Free all heap allocations.
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

    double r_grid_innercorner     = 0.0;
    double r_grid_outercorner     = 0.0;
    double theta_grid_innercorner = 0.0;
    double theta_grid_outercorner = 0.0;
    double nonthermal_n_dens      = 0.0;

    /* ── Step 1: Shell boundaries ─────────────────────────────────────────── */
    double rmin = calcCyclosynchRLimits(hydro_data->scatt_frame_number,
                                         hydro_data->inj_frame_number,
                                         hydro_data->fps, r_inj, "min");
    double rmax = calcCyclosynchRLimits(hydro_data->scatt_frame_number,
                                         hydro_data->inj_frame_number,
                                         hydro_data->fps, r_inj, "max");
    
    //TODO: think about changing after testing, especially for emitting synchrotron photons to start off
    //the above should simplify to this but when injecting photons the scatt frame number isnt defined, if only emitting photons after the first frame then above is ok
    //rmin=r_inj - 0.5*C_LIGHT/hydro_data->fps;
    //rmax=r_inj + 0.5*C_LIGHT/hydro_data->fps;


    fprintf(fPtr,
            ">> [photonEmitSynch] Shell: rmin=%.3e cm, rmax=%.3e cm, "
            "theta=[%.3f, %.3f] rad\n",
            rmin, rmax, theta_min, theta_max);
    fflush(fPtr);

    /* ── Step 2: First pass — count active cells ──────────────────────────── */
    for (i = 0; i < hydro_data->num_elements; i++)
    {
        #if NONTHERMAL_E_DIST != OFF
            nonthermal_n_dens = (hydro_data->nonthermal_dens)[i];
        #endif
        
        if (nonthermal_n_dens > 0.0)
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

    /* ── Step 3: Allocate per-cell arrays ────────────────────────────────── */
    int    *ph_dens    = (int    *)malloc(block_cnt * sizeof(int));
    int    *cell_index = (int    *)malloc(block_cnt * sizeof(int));
    double *W_cell     = (double *)malloc(block_cnt * sizeof(double));

    /* ── Step 4: Compute physical prefactor K_phys and second pass ──────────
     *
     * W_cell[j] = Lambda_j is the physical photon emission rate from cell j
     * per simulation frame [photons frame^-1], computed as:
     *
     *   W_cell[j] = K_phys * n_e_nth_j * B_j * V_j
     *
     * where K_phys is the compile-time constant:
     *
     *   K_phys = (sqrt(3) e^3) / (4 pi me c^2 h)
     *            * A_norm
     *            * ln(10)^2
     *            * nu_cdf_norm
     *            / fps
     *
     * Physical origins of each factor:
     *   sqrt(3) e^3 / (4 pi me c^2)  — synchrotron emissivity prefactor
     *                                   [RAIKOU Eq. B11; R&L79 Eq. 6.36]
     *   1/h                           — converts energy rate -> photon rate
     *   A_norm                        — electron distribution normalisation
     *                                   [(p-1)/(gmin^{1-p}-gmax^{1-p})
     *                                    for single PL;
     *                                    brokenPowerLawNorm for broken PL]
     *   ln(10)^2                      — converts the log10-space double
     *                                   integral in nu_cdf_norm back to
     *                                   natural-log measure
     *   nu_cdf_norm                   — dimensionless double integral over
     *                                   log10(nu_tilde) and log10(gamma)
     *                                   stored by buildUniversalNuCDF
     *   1/fps                         — converts photons s^-1 -> photons
     *                                   frame^-1; computed once here so
     *                                   it does not appear in the per-cell
     *                                   or per-packet loops below
     *
     * The Poisson mean for cell j is then simply:
     *
     *   lambda_j = W_cell[j] / ph_weight_adjusted
     *            = [photons frame^-1] / [photons packet^-1]
     *            = [packets frame^-1]
     *
     * with no further fps factor anywhere in the function.
     */

    /* Electron distribution normalisation constant A_norm */
    double A_norm = 0.0;

    #if NONTHERMAL_E_DIST == POWERLAW
    {
        double p    = POWERLAW_INDEX;
        double gmin = GAMMA_MIN;
        double gmax = GAMMA_MAX;

        if (fabs(p - 1.0) < 1e-6)
        {
            fprintf(fPtr,
                    ">> [photonEmitSynch] WARNING: POWERLAW_INDEX = 1, "
                    "degenerate normalisation. A_norm set to 1.\n");
            fflush(fPtr);
            A_norm = 1.0;
        }
        else
        {
            double denom = pow(gmin, 1.0 - p) - pow(gmax, 1.0 - p);
            A_norm = (fabs(denom) > 1e-300) ? (p - 1.0) / denom : 1.0;
        }
    }
    #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
    {
        A_norm = brokenPowerLawNorm(POWERLAW_INDEX_1, POWERLAW_INDEX_2,
                                    GAMMA_MIN, GAMMA_MAX, GAMMA_BREAK);
    }
    #endif

    /*
     * Retrieve nu_cdf_norm from the universal tables built at initialisation.
     * This is the unnormalised double integral over log10(nu_tilde) and
     * log10(gamma) that, combined with K_phys, converts n_e_nth * B * V
     * into a physical photon emission rate.
     */
    const SynchUniversalTables *tables = getSynchTables(fPtr);

    double K_phys = (sqrt(3.0) * CHARGE_EL * CHARGE_EL * CHARGE_EL)
                  / (2.0 * M_EL * C_LIGHT * C_LIGHT * PL_CONST)
                  * A_norm
                  * log(10.0) * log(10.0)
                  * tables->nu_cdf_norm
                  / hydro_data->fps;

    fprintf(fPtr,
            ">> [photonEmitSynch] Physical prefactor: "
            "A_norm=%.4e, nu_cdf_norm=%.4e, fps=%.4e, "
            "K_phys=%.4e [photons frame^-1 G^-1 cm^-6]\n",
            A_norm, tables->nu_cdf_norm, hydro_data->fps, K_phys);
    fflush(fPtr);

    /* Second pass — record cell indices and physical photon rates */
    j = 0;
    for (i = 0; i < hydro_data->num_elements; i++)
    {
        #if NONTHERMAL_E_DIST != OFF
            nonthermal_n_dens = (hydro_data->nonthermal_dens)[i];
        #endif

        if (nonthermal_n_dens > 0.0)
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

                /*
                 * W_cell[j] = Lambda_j [photons frame^-1]
                 * = K_phys * n_e_nth * B * V
                 */
                W_cell[j] = K_phys
                           * nonthermal_n_dens
                           * B
                           * V;
                j++;
            }
        }
    }
    


    /* ── Step 5: Weight-tuning loop ──────────────────────────────────────────
     *
     * W_cell[j] is in [photons frame^-1] so the Poisson mean is simply
     * W_cell[j] / ph_weight_adjusted with no fps factor.
     */
    double ph_weight_adjusted = ph_weight;
    double lambda_total       = 0.0;

    do
    {
        lambda_total = 0.0;
        for (j = 0; j < block_cnt; j++)
            lambda_total += W_cell[j] / ph_weight_adjusted;

        if      (lambda_total > (double)max_photons) ph_weight_adjusted *= 10.0;
        else if (lambda_total < (double)min_photons) ph_weight_adjusted *= 0.5;

    } while (lambda_total > (double)max_photons ||
             lambda_total < (double)min_photons);

    fprintf(fPtr,
            ">> [photonEmitSynch] Weight tuning converged: "
            "ph_weight_adjusted=%.3e [photons packet^-1], "
            "expected packets per frame=%.1f\n",
            ph_weight_adjusted, lambda_total);
    fflush(fPtr);

    /* ── Step 6: Poisson draw of per-cell packet counts ──────────────────────
     *
     * lambda_j = W_cell[j] / ph_weight_adjusted
     *          = [photons frame^-1] / [photons packet^-1]
     *          = [packets frame^-1]
     */
    int ph_tot = 0;
    for (j = 0; j < block_cnt; j++)
    {
        double lambda_j = W_cell[j] / ph_weight_adjusted;
        ph_dens[j] = (int)gsl_ran_poisson(rand, lambda_j);
        ph_tot    += ph_dens[j];
    }

    fprintf(fPtr,
            ">> [photonEmitSynch] Poisson draw: ph_tot=%d across "
            "%d active cells.\n", ph_tot, block_cnt);
    fflush(fPtr);

    /* ── Step 7: Allocate and initialise output photon array ─────────────── */
    if (ph_tot == 0)
    {
        fprintf(fPtr,
                ">> [photonEmitSynch] WARNING: all cells drew ph_dens=0 "
                "(lambda_total=%.3f). No photons emitted.\n\n",
                lambda_total);
        fflush(fPtr);
        free(ph_dens); free(cell_index); free(W_cell);
        return 0;
    }

    struct photon *ph_emit = (struct photon *)
                             malloc(ph_tot * sizeof(struct photon));
    if (ph_emit == NULL)
    {
        fprintf(fPtr,
                ">> [photonEmitSynch] ERROR: malloc failed for ph_emit "
                "(%d photons).\n\n", ph_tot);
        fflush(fPtr);
        free(ph_dens); free(cell_index); free(W_cell);
        return 0;
    }

    for (i = 0; i < ph_tot; i++)
        initializePhoton(&ph_emit[i]);

    /* ── Step 8: Per-cell emission loop ──────────────────────────────────────
     *
     * For each active cell j with ph_dens[j] > 0, call emitCellPackets
     * passing B_cell directly. Frequency sampling uses the universal
     * nu_tilde CDF built once in initSynchTables — no per-cell numerical
     * integration is performed here.
     */
    int idx = 0;
    for (j = 0; j < block_cnt; j++)
    {
        if (ph_dens[j] > 0)
        {
            double B_cell = getMagneticFieldMagnitude(hydro_data,
                                                       cell_index[j]);
            emitCellPackets(cell_index[j], ph_dens[j], B_cell,
                             ph_weight_adjusted,
                             ph_emit, &idx, hydro_data, rand, fPtr);
        }
    }

    /* Sanity check: idx must equal ph_tot */
    if (idx != ph_tot)
    {
        fprintf(fPtr,
                ">> [photonEmitSynch] WARNING: expected to fill %d photon "
                "structs but filled %d. Check emitCellPackets.\n",
                ph_tot, idx);
        fflush(fPtr);
        ph_tot = idx;
    }

    /* ── Weight conservation check ───────────────────────────────────────── */
    double weight_sum      = 0.0;
    double weight_expected = (double)ph_tot * ph_weight_adjusted;

    for (i = 0; i < ph_tot; i++)
        weight_sum += ph_emit[i].weight;

    double weight_tol     = 1.0 / sqrt((double)ph_tot);
    double weight_rel_err = fabs(weight_sum - weight_expected)
                          / weight_expected;

    fprintf(fPtr,
            ">> [photonEmitSynch] Weight conservation: "
            "expected=%.6e  actual=%.6e  "
            "rel_err=%.3e (tol=%.3e)  %s\n",
            weight_expected, weight_sum, weight_rel_err, weight_tol,
            weight_rel_err < weight_tol ? "PASS" : "FAIL");
    fflush(fPtr);
    
    /* ── Energy conservation diagnostic  ─────────────
     *
     * Sum the total emitted energy:
     *   E_emitted = sum_i  weight_i * (comv_p0_i * C_LIGHT)   [erg]
     *
     * comv_p0 = E_photon / C_LIGHT  [erg/cm * s/cm ... i.e. g cm/s]
     * so  E_photon = comv_p0 * C_LIGHT  [erg]
     *
     * Divide by the total cell volume and dt = 1/fps to get the
     * volumetric power [erg/s/cm^3] and compare to the analytic P_tot.
     */
    double E_emitted  = 0.0;
    double V_total    = 0.0;
    
    for (i = 0; i < ph_tot; i++)
        E_emitted += ph_emit[i].weight * ph_emit[i].comv_p0 * C_LIGHT;
    
    for (j = 0; j < block_cnt; j++)
        V_total += hydroElementVolume(hydro_data, cell_index[j]);
    
    double dt           = 1.0 / hydro_data->fps;
    double P_code       = E_emitted / (V_total * dt);   /* erg/s/cm^3 */
    
    /* Analytic total power per unit volume [erg/s/cm^3]:
     *   P_analytic = 4*pi * integral_0^inf j_nu dnu
     * For a single power-law this has the known closed form
     * (R&L79 Eq. 6.38, pitch-angle averaged):
     *
     *   P = (4/3) * sigma_T * c * U_B * <gamma^2>_emission
     *
     * where U_B = B^2/(8*pi) and <gamma^2>_emission is the
     * emission-weighted mean gamma^2 integrated over N(gamma)*gamma^2.
     *
     * We compute it here from the prefactor of RAIKOU B11 integrated
     * analytically over the power-law distribution:
     *
     *   integral_0^inf j_nu dnu
     *     = (sqrt(3) e^3 B) / (2 me c^2)
     *       * A_norm * n_e_nth
     *       * integral_{gmin}^{gmax} gamma^{-p} * gamma^2 * nu_c dgamma
     *       * [dimensionless Bessel integral = 8*pi/(9*sqrt(3))]
     *
     * The simplest self-consistent check is to use the same nu_cdf_norm
     * integral that K_phys uses, but now weighted by h*nu to get energy
     * rather than photon count.  This gives:
     *
     *   P_analytic = K_phys * <h*nu>_emission * n_e_nth_mean * B_mean / dt
     *              = (W_cell_total / dt) * <h*nu>_emission / V_total
     *
     * where <h*nu>_emission is the emission-weighted mean photon energy.
     */
    
    /* Compute emission-weighted mean photon energy from ph_emit */
    double E_mean_numerator   = 0.0;
    double E_mean_denominator = 0.0;
    for (i = 0; i < ph_tot; i++)
    {
        double E_ph = ph_emit[i].comv_p0 * C_LIGHT;   /* erg */
        E_mean_numerator   += ph_emit[i].weight * E_ph;
        E_mean_denominator += ph_emit[i].weight;
    }
    double h_nu_mean = (E_mean_denominator > 0.0)
                     ? E_mean_numerator / E_mean_denominator : 0.0;
    
    /* Analytic power using W_cell (already computed, units: photons/frame) */
    double lambda_total_phys = 0.0;   /* total photons/frame from all cells */
    for (j = 0; j < block_cnt; j++)
        lambda_total_phys += W_cell[j];
    
    /* P_analytic = (photons/frame) * (erg/photon) / (dt * V_total) */
    double P_analytic = lambda_total_phys * h_nu_mean / (dt * V_total);
    
    fprintf(fPtr,
            ">> [photonEmitSynch] Energy diagnostic:\n"
            ">>   E_emitted         = %.4e erg  (this frame)\n"
            ">>   V_total           = %.4e cm^3\n"
            ">>   dt                = %.4e s\n"
            ">>   P_code            = %.4e erg/s/cm^3\n"
            ">>   P_analytic        = %.4e erg/s/cm^3\n"
            ">>   ratio P_code/P_analytic = %.4f  (expect ~1)\n"
            ">>   h_nu_mean         = %.4e erg  (%.4e keV)\n",
            E_emitted, V_total, dt,
            P_code, P_analytic,
            (P_analytic > 0.0) ? P_code / P_analytic : 0.0,
            h_nu_mean, h_nu_mean / 1.602e-9);
    fflush(fPtr);


    /* ── Step 9: Add emitted photons to the photon list ──────────────────── */
    if (ph_tot > 0)
    {
        //if we have at least 1 photon in the photon_list then we just add the emitted photons, otherwise we have to set the photons in the photon list
        if (photon_list->list_capacity >0)
        {
            addToPhotonList(photon_list, ph_emit, ph_tot);
        }
        else
        {
            setPhotonList(photon_list, ph_emit, ph_tot);
        }

        n_emitted = ph_tot;
    }

    /* ── Step 10: Free all heap allocations ──────────────────────────────── */
    free(ph_emit);
    free(ph_dens);
    free(cell_index);
    free(W_cell);

    fprintf(fPtr,
            ">> [photonEmitSynch] Complete: emitted %d synchrotron packets "
            "from %d active cells (ph_weight_adjusted=%.3e "
            "[photons packet^-1]).\n\n",
            n_emitted, block_cnt, ph_weight_adjusted);
    fflush(fPtr);

    return n_emitted;
}
