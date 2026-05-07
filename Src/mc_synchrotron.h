/*
 * mc_synchrotron.h
 * ================
 * Header for pitch-angle averaged synchrotron photon emission with SSA.
 *
 * The physics implemented in this file follows the framework of
 * Ghisellini & Svensson (1991) [G&S91], MNRAS 252, 313, supplemented
 * by Crusius & Schlickeiser (1986) [C&S86] for the spectral kernel and
 * Kawashima et al. (2023) [K23] for the Monte Carlo weight modification.
 *
 * Equation map (G&S91 unless otherwise noted):
 * ─────────────────────────────────────────────
 *
 * Single-electron emissivity kernel R(x):
 *   R(x) = F(x) = x * integral_x^inf K_{5/3}(xi) dxi
 *   This is the standard synchrotron spectral function, see G&S91 Eq. 2
 *   and C&S86 Eq. A4. For an isotropic pitch-angle distribution the
 *   pitch-angle average reduces to F(x) (G&S91 Section 2).
 *
 * SSA absorption coefficient kappa_nu:
 *   G&S91 Eq. 12:
 *     kappa_nu = -(e^3 B) / (8 pi^2 me^2 c^2 nu^2)
 *                * integral N(gamma) * d/d(gamma)[ gamma^2 * (p_e / gamma)
 *                  * R(nu/nu_c) ] dgamma
 *   where nu_c = (3 e B gamma^2) / (4 pi me c)  [G&S91 Eq. 3]
 *
 *   After performing the gamma derivative and simplifying, G&S91 Eq. 13
 *   gives the absorption kernel as:
 *     d/d(gamma)[ gamma^2 * R(x) / (p_e * gamma) ] dx/dgamma
 *       -> 2R(x) + 2x dR/dx  [the abs_kern_arr stored in this struct]
 *
 *   The full expression for kappa_nu used in buildSynchKappaTable is
 *   G&S91 Eq. 14:
 *     kappa_nu = (sqrt(3) e^3 B) / (8 pi me c nu^2)
 *                * integral N(gamma) * [2R(x) + 2x dR/dx] dgamma
 *
 * Pitch-angle distribution:
 *   Isotropic: f(alpha) = (2/pi) sin^2(alpha)  [G&S91 Section 2, text
 *   above Eq. 2]. The factor sin(alpha) in nu_c ensures that the
 *   pitch-angle averaged emissivity integrates correctly.
 *
 * Critical frequency:
 *   nu_c = (3 e B sin(alpha) gamma^2) / (4 pi me c)  [G&S91 Eq. 3]
 *
 * Dimensionless frequency ratio:
 *   x = nu / nu_c  [G&S91 Eq. 2, argument of R(x)]
 */

#ifndef MC_SYNCHROTRON_H
#define MC_SYNCHROTRON_H

#include "mcrat.h"


/* ── Table dimension macros ───────────────────────────────────────────── */
#ifndef SYNCH_N_X
    #define SYNCH_N_X       2000   /* points in R(x) / F(x) table          */
#endif
#ifndef SYNCH_N_ALPHA
    #define SYNCH_N_ALPHA   2000   /* points in sin^2(alpha) CDF table      */
#endif
#ifndef SYNCH_N_NU
    #define SYNCH_N_NU       300   /* points in kappa_nu table per cell     */
#endif
#ifndef SYNCH_N_GL
    #define SYNCH_N_GL       200   /* Gauss-Legendre nodes for gamma integral*/
#endif
#ifndef SYNCH_N_STRATA
    #define SYNCH_N_STRATA    10   /* frequency strata for stratified sampler*/
#endif
#ifndef SYNCH_N_REF
    #define SYNCH_N_REF   500000   /* reference photons for stratum CDF      */
#endif
#ifndef SYNCH_N_CDF_BINS
    #define SYNCH_N_CDF_BINS 2000  /* bins in marginal nu CDF               */
#endif


/* ══════════════════════════════════════════════════════════════════════════
 * STRUCTURES
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * SynchUniversalTables
 * --------------------
 * All spectral tables that are independent of cell properties.
 * Allocated and initialised once via initSynchTables(); freed via
 * freeSynchTables(). Should be held at the same scope as the
 * photon_list and hydro_dataframe (e.g. at the top of the main
 * injection loop in mcrat.c).
 */
typedef struct
{
    /* x grid and R(x) = F(x) values */
    double  x_arr[SYNCH_N_X];
    double  R_arr[SYNCH_N_X];
    double  abs_kern_arr[SYNCH_N_X];   /* 2R + 2x dR/dx  [G&S91 kernel]   */

    /* Inverse CDF: uniform -> x ~ F(x) */
    double  cdf_x[SYNCH_N_X];
    gsl_spline       *inv_F_spline;
    gsl_interp_accel *inv_F_acc;

    /* Inverse CDF: uniform -> alpha ~ sin^2(alpha) */
    double  alpha_arr[SYNCH_N_ALPHA];
    double  cdf_alpha[SYNCH_N_ALPHA];
    gsl_spline       *inv_alpha_spline;
    gsl_interp_accel *inv_alpha_acc;

    /* R(x) and abs_kern(x) forward interpolators */
    gsl_spline       *R_spline;
    gsl_interp_accel *R_acc;
    gsl_spline       *abs_kern_spline;
    gsl_interp_accel *abs_kern_acc;

} SynchUniversalTables;

/*
 * SynchKappaTable
 * ---------------
 * Per-cell SSA opacity table: log10(kappa_nu) vs log10(nu).
 * Built for the representative cell (B_ref, ne_ref) and rescaled
 * linearly for other cells since kappa_nu ∝ n_e * B.
 */
typedef struct
{
    int     n_nu;
    double *log_nu;      /* log10(nu) [n_nu]          */
    double *log_kappa;   /* log10(kappa_nu) [n_nu]    */
    double  B_ref;       /* B used to build this table [G]    */
    double  ne_ref;      /* n_e used to build this table [cm^-3] */
    gsl_spline       *spline;
    gsl_interp_accel *acc;
} SynchKappaTable;

/*
 * SynchStratifiedParams
 * ---------------------
 * Parameters for stratified frequency sampling built from a reference
 * sample of the natural emission distribution.
 */
typedef struct
{
    int    n_strata;
    double strata_edges[SYNCH_N_STRATA + 1];
    double stratum_probs[SYNCH_N_STRATA];    /* P(nu in stratum k)          */
    /* Marginal CDF of log10(nu) */
    double cdf_log_nu_edges[SYNCH_N_CDF_BINS + 1];
    double cdf_log_nu_vals[SYNCH_N_CDF_BINS + 1];
    gsl_spline       *inv_nu_cdf_spline;
    gsl_interp_accel *inv_nu_cdf_acc;
} SynchStratifiedParams;

/* ══════════════════════════════════════════════════════════════════════════
 * FUNCTION PROTOTYPES
 * ══════════════════════════════════════════════════════════════════════════ */

/* Initialise / free universal spectral tables */
void initSynchTables(SynchUniversalTables *tables, FILE *fPtr);
void freeSynchTables(SynchUniversalTables *tables);

/* Build / free per-cell kappa_nu table */
SynchKappaTable *buildSynchKappaTable(double B,
                                       const SynchUniversalTables *tables,
                                       FILE *fPtr);
void             freeSynchKappaTable(SynchKappaTable *kt);

/* Build / free stratified frequency sampler */
void buildSynchStratifiedParams(SynchStratifiedParams *sp,
                                 const SynchUniversalTables *tables,
                                 struct hydro_dataframe *hydro_data,
                                 double nu_min, double nu_max,
                                 gsl_rng *rand, FILE *fPtr);
void freeSynchStratifiedParams(SynchStratifiedParams *sp);

/*
 * photonEmitSynch
 * ---------------
 * Main emission function. Signature mirrors photonEmitCyclosynch so it
 * can be called from the same site in mcrat.c.
 *
 * Parameters
 * ----------
 * photon_list  : live photon list to append new photons into
 * r_inj        : injection radius [cm]
 * ph_weight    : suggested photon weight
 * maximum_photons : upper bound on photons to emit this call
 * theta_min/max: jet opening angle range [rad]
 * hydro_data   : fluid grid
 * tables       : pre-built universal spectral tables
 * sp           : pre-built stratified sampler parameters
 * rand         : GSL RNG
 * fPtr         : log file pointer
 *
 * Returns number of photons added to photon_list.
 */
int photonEmitSynch(struct photonList       *photon_list,
                    double                   r_inj,
                    double                   ph_weight,
                    int                      maximum_photons,
                    double                   theta_min,
                    double                   theta_max,
                    struct hydro_dataframe  *hydro_data,
                    const SynchUniversalTables *tables,
                    const SynchStratifiedParams *sp,
                    gsl_rng                 *rand,
                    FILE                    *fPtr);

#endif /* MC_SYNCHROTRON_H */
