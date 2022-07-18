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
