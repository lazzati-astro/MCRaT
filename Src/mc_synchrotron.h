/*
 * mc_synchrotron.h
 * ================
 * Pitch-angle averaged synchrotron photon emission and synchrotron
 * self-absorption (SSA) for non-thermal power-law and broken power-law
 * electron distributions.
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * PHYSICS OVERVIEW
 * ─────────────────────────────────────────────────────────────────────────────
 *
 * EMISSION
 * --------
 * The pitch-angle averaged synchrotron emissivity integrated over an
 * isotropic pitch-angle distribution for a single power-law electron
 * distribution is (RAIKOU Eq. B11; G&S91 Eq. 2.16):
 *
 *   j_{nu_f}^(f) = (p-1) n_{e,nth} e^2 nu_c
 *                  / (2 sqrt(3) c (gamma_min^{1-p} - gamma_max^{1-p}))
 *                  * (nu_f/nu_c)^{-(p-1)/2}
 *                  * [G(x_max) - G(x_min)]
 *
 * where
 *   G(x)  = integral_x^inf z^{(p-3)/2} F(z) dz   [RAIKOU Eq. B12]
 *   F(z)  = z * integral_z^inf K_{5/3}(xi) dxi    [RAIKOU Eq. B13]
 *   nu_c  = 3 e B / (4 pi me c)  (pitch-angle averaged)  [RAIKOU Eq. B7]
 *   x     = nu_f / (gamma_e^2 nu_c)
 *
 * ABSORPTION
 * ----------
 * The SSA absorption coefficient in the fluid rest frame for a single
 * power-law electron distribution is (RAIKOU Eq. C2):
 *
 *   alpha_{nu_f}^(f) = (p-1)(p+2) n_{e,nth} e^2 nu_c
 *                      / (4 sqrt(3) me c (gamma_min^{1-p} - gamma_max^{1-p}))
 *                      * (nu_f/nu_c)^{-(p+4)/2}
 *                      * [G_a(x_max) - G_a(x_min)]
 *
 * where the universal absorption spectral integral is (RAIKOU Eq. C3):
 *
 *   G_a(x) = integral_x^inf z^{(p-2)/2} F(z) dz
 *
 * For a broken power-law the equivalent is RAIKOU Eq. C4, summing two
 * segment contributions with indices p1 and p2 respectively.
 *
 * G_a(x) depends on p and x = nu/(nu_c * gamma^2) but NOT on B or nu_f
 * individually. A single 1D table of G_a(x) therefore serves all cells,
 * with B and n_e entering analytically through the RAIKOU Eq. C2 prefactor
 * at runtime. This is the RAIKOU analytic approach described in Appendix C.
 *
 * WEIGHT MODIFICATION
 * -------------------
 * The photon weight is attenuated by SSA along each path step dl following
 * RAIKOU Eq. 40:
 *
 *   w_new = w_old * exp(-Delta tau_{nu_f}^(a))
 *
 * where the optical depth increment is (RAIKOU Eq. 31):
 *
 *   Delta tau_{nu_f}^(a) = (nu_f / nu_z) * alpha_{nu_f}^(f) * dl
 *
 * ─────────────────────────────────────────────────────────────────────────────
 * REFERENCES
 * ─────────────────────────────────────────────────────────────────────────────
 *   RAIKOU : Kawashima, Ohsuga & Takahashi (2023), ApJ 949 101
 *   G&S91  : Ghisellini & Svensson (1991), MNRAS 252, 313
 *   R&L79  : Rybicki & Lightman (1979), Radiative Processes in Astrophysics
 */

#ifndef MC_SYNCHROTRON_H
#define MC_SYNCHROTRON_H

/* ── Table dimensions ──────────────────────────────────────────────────────── */

/*
 * Number of log-spaced points in the universal x-grid used for both F(x)
 * [RAIKOU Eq. B13] and G_a(x) [RAIKOU Eq. C3]. 2000 points gives
 * interpolation errors below 1e-4 across the dynamic range of interest.
 */
#define SYNCH_N_X        2000

/*
 * Bounds of the x-grid for F(x) and G_a(x).
 * F(x) decays exponentially for x > 50; below 1e-5 the power-law asymptote
 * is well represented by the spline boundary condition.
 */
#define SYNCH_X_MIN      1e-5
#define SYNCH_X_MAX      50.0

/*
 * Number of log-frequency strata for stratified importance sampling.
 * One stratum per decade of frequency, spanning nu_c(gamma_min) to
 * nu_c(gamma_max).
 */
#define SYNCH_N_STRATA   10

/*
 * Number of photons in the reference sample used to measure stratum
 * probabilities in buildSynchCellStrata. 500000 gives estimates
 * accurate to ~0.1%.
 */
#define SYNCH_N_REF      500000

/*
 * Number of fine-grained bins in the marginal CDF of log10(nu_f) stored
 * inside SynchStratifiedParams.
 */
#define SYNCH_N_CDF_BINS 1000

/*
 * Minimum stratum probability below which the stratum is skipped during
 * photon emission.
 */
#define SYNCH_P_MIN      1e-7


/* ── Universal spectral tables ─────────────────────────────────────────────── */

/*
 * SynchUniversalTables
 * --------------------
 * All spectral functions that depend only on the dimensionless frequency
 * ratio x = nu / (gamma^2 nu_c), built once at initialisation and shared
 * across all cells and photons.
 *
 * F(x)    : synchrotron spectral function  [RAIKOU Eq. B13]
 * G_a(x)  : absorption spectral integral  [RAIKOU Eq. C3]
 *            For a single power law one G_a table is built with POWERLAW_INDEX.
 *            For a broken power law two tables are built: one with
 *            POWERLAW_INDEX_1 (Ga_spline_p1) and one with POWERLAW_INDEX_2
 *            (Ga_spline_p2), matching the two segment contributions in
 *            RAIKOU Eq. C4.
 *
 * Inverse CDF of x ~ F(x)*x d(log x): used to sample x for emission
 *   [G&S91 Sec. 2; R&L79 Eq. 6.36]
 *
 * Inverse CDF of alpha ~ sin^2(alpha): used to sample pitch angle
 *   [G&S91 Sec. 2; R&L79 Sec. 6.2]
 */
typedef struct
{
    /* Shared x-grid */
    double  *x_arr;              /* [SYNCH_N_X] log-spaced x values          */

    /* F(x)  [RAIKOU Eq. B13] */
    double  *F_arr;              /* [SYNCH_N_X] F(x) values                  */
    gsl_spline       *F_spline;
    gsl_interp_accel *F_acc;

    /* G_a(x; p)  [RAIKOU Eq. C3] */
    double  *Ga_arr;             /* [SYNCH_N_X]  single PL or alias for p1   */
    double  *Ga_arr_p1;          /* [SYNCH_N_X]  broken PL low segment       */
    double  *Ga_arr_p2;          /* [SYNCH_N_X]  broken PL high segment      */
    gsl_spline       *Ga_spline;
    gsl_spline       *Ga_spline_p1;
    gsl_spline       *Ga_spline_p2;
    gsl_interp_accel *Ga_acc;
    gsl_interp_accel *Ga_acc_p1;
    gsl_interp_accel *Ga_acc_p2;

    /* Inverse CDF of x ~ F(x)*x  (emission x sampler)  [G&S91 Sec. 2] */
    double  *inv_x_cdf_u;        /* [SYNCH_N_X] CDF values u in [0,1]        */
    double  *inv_x_cdf_logx;     /* [SYNCH_N_X] corresponding log10(x)       */
    gsl_spline       *inv_x_cdf_spline;
    gsl_interp_accel *inv_x_cdf_acc;

    /* Inverse CDF of alpha ~ sin^2(alpha) (pitch-angle sampler)         */
    double  *inv_alpha_cdf_u;    /* [SYNCH_N_X] CDF values u in [0,1]        */
    double  *inv_alpha_cdf_alpha;/* [SYNCH_N_X] corresponding alpha in [0,pi]*/
    gsl_spline       *inv_alpha_cdf_spline;
    gsl_interp_accel *inv_alpha_cdf_acc;

} SynchUniversalTables;


/* ── Stratified frequency sampler parameters ───────────────────────────────── */

/*
 * SynchCellStrata
 * ---------------
 * Precomputed conditional frequency sampling parameters for a single fluid
 * cell, built from that cell's own magnetic field B_cell.
 *
 * For a cell with critical frequency nu_c = 3 e B_cell / (4 pi me c)
 * [RAIKOU Eq. B7], the natural synchrotron emission spectrum spans
 * [nu_c * gamma_min^2, nu_c * gamma_max^2]. This range is divided into
 * SYNCH_N_STRATA equal log-frequency intervals (strata). The probability
 * p_k that a naturally emitted photon falls in stratum k is measured from
 * a reference sample of SYNCH_N_REF draws using synchNaturalNu at B_cell.
 *
 * During emission, n_cell photons are divided equally across active strata.
 * Each photon in stratum k carries the importance weight correction
 *   w_k = p_k * SYNCH_N_STRATA
 * so the weighted spectrum recovers RAIKOU Eq. B11 for this cell's B_cell.
 *
 * The inverse marginal CDF spline maps u in [0,1] -> log10(nu_f) and is
 * used to sample nu_f within a stratum without rejection.
 */
typedef struct
{
    double  B_cell;                               /* field this was built for [G] */
    double  strata_edges[SYNCH_N_STRATA + 1];     /* nu boundaries [Hz]          */
    double  stratum_probs[SYNCH_N_STRATA];        /* p_k from reference sample   */

    /* Fine-grained marginal CDF of log10(nu_f) */
    double  cdf_log_nu_edges[SYNCH_N_CDF_BINS + 1];
    double  cdf_log_nu_vals [SYNCH_N_CDF_BINS + 1];

    /* Inverse CDF spline: u -> log10(nu_f) */
    gsl_spline       *inv_nu_cdf_spline;
    gsl_interp_accel *inv_nu_cdf_acc;

} SynchCellStrata;


/* ── Function prototypes ───────────────────────────────────────────────────── */

/* Universal table lifecycle */
void   initSynchTables (SynchUniversalTables *tables, FILE *fPtr);
void   freeSynchTables (SynchUniversalTables *tables);

/* Emission MC samplers */
double synchSampleX            (const SynchUniversalTables *tables, double u);
double synchSampleAlpha        (const SynchUniversalTables *tables, double u);
double synchSampleGammaEmission(gsl_rng *rand);

/* Absorption coefficient  [RAIKOU Appendix C] */
double synchAlphaNu(double nu_f,
                    double B,
                    double n_e_nth,
                    const SynchUniversalTables *tables,
                    FILE *fPtr);

/* SSA weight modification  [RAIKOU Eqs. 31, 40] */
void applySynchSSAWeightModification(struct photon              *ph,
                                      double                      dl,
                                      double                      B_cell,
                                      double                      n_e_nth_cell,
                                      const SynchUniversalTables *tables,
                                      FILE                       *fPtr);

/* Stratified sampler lifecycle */
void buildSynchCellStrata  (SynchCellStrata            *cs,
                             double                      B_cell,
                             const SynchUniversalTables *tables,
                             gsl_rng                    *rand,
                             FILE                       *fPtr);
void freeSynchCellStrata   (SynchCellStrata *cs);

/* Main emission function */
int photonEmitSynch(struct photonList          *photon_list,
                    double                      r_inj,
                    double                      ph_weight,
                    int                         min_photons,
                    int                         max_photons,
                    double                      theta_min,
                    double                      theta_max,
                    struct hydro_dataframe     *hydro_data,
                    const SynchUniversalTables *tables,
                    gsl_rng                    *rand,
                    FILE                       *fPtr);

#endif /* MC_SYNCHROTRON_H */
