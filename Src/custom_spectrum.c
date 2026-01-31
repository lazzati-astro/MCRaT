//
//  custom_spectrum.c
//  
//
//  Created by Tyler Parsotan on 7/14/22.
//
//  This capability is currently untested
//  MCRAT calls the custom_photon_sampler function which allows a user to overwrite the photons that have been initalized in the simulation.
//  Currently, the code allows the user to only modify the photons comoving or lab 4 momentum, or the stokes parameters. If the user wants to modify the energy of the injected photon (while maintaining the presampled angular directions) they can just set the p0 or comv_p0 values of the 4 momentum.
//  If either the lab or the fluid photon 4 momentum is specified the code will properly calculate the other frame's 4 momentum doing the appropriate lorentz transform.
//  The code expects certain normalizations/units for the 4 momentum (which can be found in the MCRaT PDF documentation). As for the stokes, the code expects s0 to be 1 and s1=q=Q/I, etc. It also checks that s0^2 >= s1^2 + s2^2 + s3^2

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

struct photon custom_photon_sampler(struct hydro_dataframe *hydro_data, int hydro_index, gsl_rng * rand, FILE *fPtr)
{
    struct photon custom_photon=createPhoton();
    
    double bb_temp=1e-8 *(GSL_CONST_CGSM_MASS_ELECTRON*GSL_CONST_CGSM_SPEED_OF_LIGHT*GSL_CONST_CGSM_SPEED_OF_LIGHT)/GSL_CONST_CGSM_BOLTZMANN;
    double fr_dum;
    
    double fr_max=(3.31e10)*bb_temp;//max frequency of bb photon density spectrum
    double bb_norm=( pow((fr_max),2.0))/gsl_expm1(PL_CONST*fr_max/(K_B*bb_temp));  //find value of bb at fr_max
    double y_dum=1; //initalize loop
    double yfr_dum=0;
    while (y_dum>yfr_dum)
    {
        fr_dum=gsl_rng_uniform_pos(rand)*6.3e11*bb_temp; //in Hz
        //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
        y_dum=gsl_rng_uniform_pos(rand);
        
        yfr_dum=((1.0/bb_norm)* pow((fr_dum),2.0))/gsl_expm1(PL_CONST*fr_dum/(K_B*bb_temp)); //(exp(PL_CONST*fr_dum/(K_B*bb_temp))-1); //curve is normalized to vaue of bb @ max frequency
    }
    
    //just set the energy since the code will recalculate the 4 mometum, assuming anisotropic photon angle distribution. The code will also recalculate the lab frame photon 4 momentum by doing the proper lorentz boost
    custom_photon.comv_p0=fr_dum*GSL_CONST_CGSM_PLANCKS_CONSTANT_H/GSL_CONST_CGSM_SPEED_OF_LIGHT;
    

    return custom_photon;
    
}
