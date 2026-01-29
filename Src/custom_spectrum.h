//
//  custom_spectrum.h
//  
//
//  Created by Tyler Parsotan on 7/14/22.
//

double custom_spectrum(double frequency);

struct photon custom_photon_sampler(struct hydro_dataframe *hydro_data, int hydro_index, gsl_rng * rand, FILE *fPtr);
