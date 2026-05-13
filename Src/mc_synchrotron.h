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
 * FREQUENCY SAMPLING — UNIVERSAL CDF
 * ------------------------------------
 * The photon-number emissivity per unit log-frequency is:
 *
 *   dN/d(ln nu) propto integral gamma^{1-p} F(nu/(gamma^2 nu_c)) d(ln gamma)
 *
 * Substituting the dimensionless frequency nu_tilde = nu / nu_c:
 *
 *   dN/d(ln nu_tilde) propto integral gamma^{1-p} F(nu_tilde/gamma^2) d(ln gamma)
 *
 * The magnetic field B has disappeared completely. The CDF shape depends
 * only on p, gamma_min, gamma_max — all compile-time constants. A single
 * universal CDF in nu_tilde is therefore built once at initialisation inside
 * initSynchTables. At emission time the physical frequency is recovered as:
 *
 *   nu_f = nu_tilde * nu_c(B_cell)   where nu_c = 3 e B / (4 pi me c)
 *
 * This eliminates the per-cell numerical integration that was previously
 * required, reducing the cost from O(N_cells * N_nu * N_gamma) per emission
 * epoch to O(1) per photon.
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
 * at runtime.
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

/* ─────────────────────────────────────────────────────────────────────────────
 * TABLE DIMENSIONS
 * ─────────────────────────────────────────────────────────────────────────────
 */

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
 * Number of log-spaced points in the universal dimensionless-frequency
 * (nu_tilde = nu/nu_c) grid used to build the photon-number CDF.
 * 2000 points gives sub-percent interpolation error across the full dynamic
 * range [SYNCH_X_MIN * GAMMA_MIN^2, SYNCH_X_MAX * GAMMA_MAX^2].
 */
#define SYNCH_N_NU        2000

/*
 * Number of log-gamma integration points used in the trapezoid rule over
 * the electron distribution when building the universal nu_tilde CDF.
 * 300 points gives < 0.1% error on a power-law integrand.
 */
#define SYNCH_N_GAMMA_INT  300

/*
 * Number of equal log-nu_tilde strata used for stratified frequency
 * sampling. Each stratum receives an equal number of photon packets;
 * the importance weight w_k = w_0 * p_k * K corrects for the forced
 * equal allocation. ~1 stratum per decade of frequency is recommended.
 */
#define SYNCH_N_STRATA      20


/* ─────────────────────────────────────────────────────────────────────────────
 * UNIVERSAL SPECTRAL TABLES
 * ─────────────────────────────────────────────────────────────────────────────
 *
 * SynchUniversalTables
 * --------------------
 * All spectral quantities that depend only on compile-time constants
 * (p, gamma_min, gamma_max) or the dimensionless ratio x = nu/(gamma^2 nu_c).
 * Built once by initSynchTables at startup and shared read-only across all
 * cells, photons, and OpenMP threads.
 *
 * Thread safety
 * -------------
 * The splines themselves are read-only after initSynchTables returns.
 * The mutable gsl_interp_accel search caches are NOT stored here — they
 * live in a per-thread SynchThreadAccel array inside mc_synchrotron.c.
 *
 * Members
 * -------
 * x_arr / F_arr / F_spline
 *   Log-spaced x-grid and the synchrotron spectral function F(x)
 *   [RAIKOU Eq. B13] evaluated on it.
 *
 * Ga_arr / Ga_arr_p1 / Ga_arr_p2
 *   G_a(x; p) [RAIKOU Eq. C3] evaluated on x_arr.
 *   For POWERLAW:      Ga_arr holds G_a at POWERLAW_INDEX.
 *   For BROKENPOWERLAW: Ga_arr_p1 holds G_a at POWERLAW_INDEX_1,
 *                       Ga_arr_p2 holds G_a at POWERLAW_INDEX_2,
 *                       Ga_arr is an alias for Ga_arr_p1.
 *   Ga_spline / Ga_spline_p1 / Ga_spline_p2 are the corresponding splines.
 *
 * inv_x_cdf_u / inv_x_cdf_logx / inv_x_cdf_spline
 *   Inverse CDF of x ~ F(x)*x d(log x), stored in log-complementary space
 *   v = -log10(1-u) for numerical stability near u -> 1. Used to sample
 *   the single-electron emission x directly [G&S91 Sec. 2; R&L79 Eq. 6.36].
 *
 * inv_alpha_cdf_u / inv_alpha_cdf_alpha / inv_alpha_cdf_spline
 *   Inverse CDF of the isotropic pitch-angle distribution
 *   f(alpha) = (2/pi) sin^2(alpha) [G&S91 Sec. 2; R&L79 Sec. 6.2].
 *
 * nu_cdf_u / nu_cdf_log_nu_tilde / n_nu_cdf
 *   Universal photon-number CDF in nu_tilde = nu/nu_c space. Stores only
 *   strictly-increasing (u, log10(nu_tilde)) pairs. At emission time the
 *   physical frequency is nu_f = 10^(sampled log10(nu_tilde)) * nu_c(B_cell).
 *
 * strata_log_nu_tilde_edges / strata_cdf_lo / strata_cdf_hi / strata_p_k
 *   Stratified sampling parameters. The nu_tilde range is divided into
 *   SYNCH_N_STRATA equal log intervals. strata_p_k[k] is the fraction of
 *   the total photon-number spectrum that falls in stratum k, used as the
 *   importance weight correction w_k = w_0 * p_k * SYNCH_N_STRATA.
 */
typedef struct
{
    /* ── x-grid ────────────────────────────────────────────────────────── */
    double  *x_arr;                  /* [SYNCH_N_X] log-spaced x values     */

    /* ── F(x)  [RAIKOU Eq. B13] ────────────────────────────────────────── */
    double  *F_arr;                  /* [SYNCH_N_X] F(x) values             */
    gsl_spline *F_spline;

    /* ── G_a(x; p)  [RAIKOU Eq. C3] ────────────────────────────────────── */
    double  *Ga_arr;                 /* single PL, or alias for Ga_arr_p1   */
    double  *Ga_arr_p1;              /* broken PL low-energy segment        */
    double  *Ga_arr_p2;              /* broken PL high-energy segment       */
    gsl_spline *Ga_spline;
    gsl_spline *Ga_spline_p1;
    gsl_spline *Ga_spline_p2;

    /* ── Inverse CDF of x ~ F(x)*x  [G&S91 Sec. 2] ─────────────────────── */
    double  *inv_x_cdf_u;            /* v = -log10(1-u) values              */
    double  *inv_x_cdf_logx;         /* corresponding log10(x)              */
    gsl_spline *inv_x_cdf_spline;

    /* ── Inverse CDF of alpha ~ sin^2(alpha)  [G&S91 Sec. 2] ───────────── */
    double  *inv_alpha_cdf_u;        /* u in [0,1]                          */
    double  *inv_alpha_cdf_alpha;    /* corresponding alpha in [0, pi]      */
    gsl_spline *inv_alpha_cdf_spline;

    /* ── Universal photon-number CDF in nu_tilde space ──────────────────── */
    double   nu_cdf_u           [SYNCH_N_NU]; /* u in [0,1]                 */
    double   nu_cdf_log_nu_tilde[SYNCH_N_NU]; /* log10(nu_tilde)            */
    int      n_nu_cdf;                        /* number of valid CDF points */
    
    /*
     * nu_cdf_norm : unnormalised integral of the pdf before normalisation,
     *
     *   nu_cdf_norm = integral_{log_nt_min}^{log_nt_max}
     *                   integral_{log_g_min}^{log_g_max}
     *                     gamma^{1-p} F(nu_tilde/gamma^2)
     *                   d(log10 gamma) d(log10 nu_tilde)
     *
     * Combined with the physical prefactor K_phys (see photonEmitSynch),
     * this converts W_cell = n_e_nth * B * V into the actual photon
     * emission rate Lambda [photons s^-1] for each cell.
     */
    double   nu_cdf_norm;

    /* ── Stratified sampling boundaries ─────────────────────────────────── */
    double   strata_log_nu_tilde_edges[SYNCH_N_STRATA + 1]; /* log10(nu_tilde) edges  */
    double   strata_cdf_lo            [SYNCH_N_STRATA];     /* CDF at lower edge      */
    double   strata_cdf_hi            [SYNCH_N_STRATA];     /* CDF at upper edge      */
    double   strata_p_k               [SYNCH_N_STRATA];     /* = cdf_hi - cdf_lo      */

} SynchUniversalTables;


/* ─────────────────────────────────────────────────────────────────────────────
 * FUNCTION PROTOTYPES
 * ─────────────────────────────────────────────────────────────────────────────
 */

/* ── Table lifecycle ────────────────────────────────────────────────────────
 *
 * initSynchTables  — build all universal tables; call once at startup before
 *                    any photon emission or SSA evaluation.
 * freeSynchTables  — release all heap memory; call once at end of simulation.
 */
void initSynchTables(FILE *fPtr);
void freeSynchTables(FILE *fPtr);

/* ── Frequency sampler ──────────────────────────────────────────────────────
 *
 * sampleSynchFrequency
 *   Draw a photon frequency from stratum k of the universal nu_tilde CDF
 *   and scale to physical Hz via nu_c(B_cell) = 3 e B / (4 pi me c).
 *
 *   Parameters
 *     k       : stratum index in [0, SYNCH_N_STRATA)
 *     B_cell  : magnetic field magnitude in the emitting cell [G]
 *     rand    : GSL RNG
 *
 *   Returns nu_f [Hz] in the fluid rest frame.
 */
double sampleSynchFrequency(int k, double B_cell, gsl_rng *rand);

/* ── Absorption coefficient  [RAIKOU Appendix C] ───────────────────────────
 *
 * synchAlphaNu
 *   Compute the SSA absorption coefficient alpha_{nu_f}^(f) [cm^-1] in the
 *   fluid rest frame at comoving frequency nu_f [Hz].
 *
 *   For POWERLAW      uses RAIKOU Eq. C2.
 *   For BROKENPOWERLAW uses RAIKOU Eq. C4.
 */
double synchAlphaNu(double nu_f,
                    double B,
                    double n_e_nth,
                    FILE  *fPtr);

/* ── SSA weight modification  [RAIKOU Eqs. 31, 40] ─────────────────────────
 *
 * calculateOpticalDepthSSA
 *   Store alpha_{nu_f}^(f) in ph->abs_optical_depth for later use by
 *   applyabsorption during photon transport. Called once per cell crossing.
 *
 * applyabsorption
 *   Attenuate ph->weight by exp(-abs_optical_depth * dl) over a lab-frame
 *   step of length dl [cm]. Called at each transport sub-step.
 *
 * applySynchSSAWeightModification
 *   Combined convenience function: compute alpha and apply the full
 *   RAIKOU Eq. 31 frame correction (nu_f/nu_z) * alpha * dl in one call.
 *   Used at photon birth to pre-attenuate the weight for the birth-cell
 *   path length.
 */
void calculateOpticalDepthSSA(struct photon          *ph,
                               struct hydro_dataframe *hydro_data,
                               FILE                   *fPtr);

void applyabsorption(struct photon *ph, double dl);

void applySynchSSAWeightModification(struct photon              *ph,
                                      double                      dl,
                                      double                      B_cell,
                                      double                      n_e_nth_cell,
                                      FILE                       *fPtr);

/* ── Main emission function ─────────────────────────────────────────────────
 *
 * photonEmitSynch
 *   Emit synchrotron photon packets from all active cells in the injection
 *   shell at radius r_inj. Returns the number of packets emitted.
 *
 *   Active cells satisfy:
 *     (a) nonthermal_dens > 0
 *     (b) cell bounding box overlaps [rmin, rmax] in radius
 *     (c) cell bounding box overlaps [theta_min, theta_max] in polar angle
 *
 *   Packet counts are drawn from Poisson distributions with means
 *   proportional to n_e_nth * B^2 * V, tuned so that the total lies in
 *   [min_photons, max_photons]. Frequencies are sampled from the universal
 *   nu_tilde CDF using stratified sampling with SYNCH_N_STRATA strata.
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
                    FILE                       *fPtr);

#endif /* MC_SYNCHROTRON_H */
