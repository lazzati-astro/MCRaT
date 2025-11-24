//
// Created by Tyler Parsotan on 11/8/25.
//

void calculateOpticalDepth(struct photon *ph, struct hydro_dataframe *hydro_data, gsl_rng *rand, FILE *fPtr);

void getCrossSection(double photon_comv_e, double fluid_temp, double *cross_section, gsl_rng *rand, FILE *fPtr);

double getThermalCrossSection(double photon_comv_e, double fluid_temp, gsl_rng *rand, FILE *fPtr);

#if NONTHERMAL_E_DIST != OFF
    double getNonThermalCrossSection(double photon_comv_e, double *subgroup_interpolated_results, gsl_rng *rand, FILE *fPtr);
#endif