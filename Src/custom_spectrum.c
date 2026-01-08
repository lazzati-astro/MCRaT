//
//  custom_spectrum.c
//  
//
//  Created by Tyler Parsotan on 7/14/22.
//

#include "mcrat.h"
#include <gsl/gsl_sf_exp.h>

double custom_spectrum(double frequency)
{
    double std=0.2*1.60218e-9/PL_CONST, mean=6.4*1.60218e-9/PL_CONST; //variance and mean in keV converted to be in Hz
    double val=0;
    
    if ((frequency>mean-6*std) && (frequency<mean+6*std))
    {
        val=gsl_sf_exp(-(frequency-mean)*(frequency-mean)/(2*std*std));//dont need (1/(std*sqrt(2*M_PI))) factor since its already normalized
    }
    else
    {
        val=0;
    }
    
    return val;
}

double custom_photon_sampler(struct hydro_dataframe *hydro_data, int hydro_index, gsl_rng * rand, FILE *fPtr)
{
    struct photon custom_photon;
    
    double bb_temp=1e-8 *(GSL_CONST_CGSM_MASS_ELECTRON*GSL_CONST_CGSM_SPEED_OF_LIGHT*GSL_CONST_CGSM_SPEED_OF_LIGHT)/GSL_CONST_CGSM_BOLTZMANN;
    double test=0;
    double test_rand1=gsl_rng_uniform_pos(rand);
    double test_rand2=gsl_rng_uniform_pos(rand);
    double test_rand3=gsl_rng_uniform_pos(rand);
    double test_rand4=gsl_rng_uniform_pos(rand);
    double test_rand5=gsl_rng_uniform_pos(rand);
    double test_cnt=0;
    double fr_dum;
    
    while (test<M_PI*M_PI*M_PI*M_PI*test_rand1/90.0)
    {
        test_cnt+=1;
        test+=1/(test_cnt*test_cnt*test_cnt*test_cnt);
    }
    fr_dum=-log(test_rand2*test_rand3*test_rand4*test_rand5)/test_cnt;
    fr_dum*=K_B*bb_temp/PL_CONST;

    return fr_dum;
    
}
