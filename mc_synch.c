/*
This file is for the different functions for emitting and absorbing synchrotron photons
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <glob.h>
#include <unistd.h>
#include <dirent.h>
#include "hdf5.h"
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include "mclib.h"
#include <omp.h>
#include "mpi.h"

//#DEFINE CRITICAL_B FINE_STRUCT*sqrt(M_EL*C_LIGHT*C_LIGHT/pow(R_EL, 3))

double calcCyclotronFreq(double magnetic_field)
{
    //B has to be in gauss
    return CHARGE_EL*magnetic_field/(2*M_PI*M_EL*C_LIGHT);
}

double calcDimlessTheta(double temp)
{
    //temp has to be in kelvin
    return K_B*temp/(M_EL*C_LIGHT*C_LIGHT);
}

double calcB(double el_dens, double temp, double epsilon_b)
{
    //calc the B field from assuming its some fraction of the matter energy density
    return sqrt(8*M_PI*3*el_dens*K_B*temp/2);
}

double n_el_MJ(double el_dens, double dimlesstheta, double gamma)
{
    //function to calulate the number density of electrons using the maxwell juttner distribution 
    return el_dens*gamma*sqrt(gamma*gamma-1)*exp(-gamma/dimlesstheta)/(dimlesstheta*gsl_sf_bessel_Kn(2, 1.0/dimlesstheta));
}

double n_el_MB(double el_dens, double dimlesstheta, double gamma)
{
    //function to calc the number density of electrons at a given dimensionless temp and lorentz factor with the maxwell boltzmann dist
    double temp=dimlesstheta*(M_EL*C_LIGHT*C_LIGHT)/K_B;
    double v=C_LIGHT*sqrt(1-(1/pow(gamma, 2)));
    
    return el_dens*4*M_PI*pow(M_EL/(2*M_PI*K_B*temp) , 3/2)*(v*C_LIGHT*C_LIGHT/(pow(gamma, 3)))*exp((-M_EL*pow(v, 2))/(2*K_B*temp));
}

//These functions are to calculate the emissivity from Wardzinski+ 2000
double Z(double nu, double nu_c, double gamma )
{
    return pow(sqrt(pow(gamma,2)-1)*np.exp(1/gamma)/(1+gamma) ,2*nu*gamma/nu_c);
}

double Z_sec_der(double nu, double nu_c, double gamma)
{
    //calculated from mathematica and plugged in theta (from paper=pi/2)
    return nu*(-2*pow(gamma,3)*(1+gamma) + 4*pow(gamma,4)*(1+gamma-pow(gamma,2)-pow(gamma,3))*log(sqrt(pow(gamma,2)-1)*exp(1/gamma)/(1+gamma) ))/(nu_c*pow(gamma,5)*(1+gamma));
}

double chi(double dimlesstheta, double gamma)
{
    double val=0;
    
    if (dimlesstheta<=0.08)
    {
        val=sqrt(2*dimlesstheta*(pow(gamma,2)-1)/(gamma*(3*pow(gamma,2)-1)));
    }
    else
    {
        val=sqrt(2*dimlesstheta/(3*gamma));
    }
    
    return val;
}

double gamma0(double nu, double nu_c, double dimlesstheta)
{
    double val=0;

    if (dimlesstheta<=0.08)
    {
        val=sqrt(pow(1+(2*nu*dimlesstheta/nu_c)*(1+(9*nu*dimlesstheta/(2*nu_c))), (-1.0/3.0) ));
    }
    else
    {
        val=sqrt(pow((1+(4*nu*dimlesstheta/(3*nu_c))), (2.0/3.0)) );
    }
    
    return val;
}

double jnu(double nu, double nu_c, double dimlesstheta, double el_dens)
{
    double dimlesstheta_ref=calcDimlessTheta(1e7);
    double gamma=gamma0(nu, nu_c, dimlesstheta);
    double val=0;
    
    if (dimlesstheta<dimlesstheta_ref)
    {
        val=(pow(M_PI,(3.0/2.0))*pow(CHARGE_EL, 2)/(pow(2,(3.0/2.0))*C_LIGHT))*sqrt(nu*nu_c)*n_el_MB(el_dens, dimlesstheta, gamma)* Z(nu, nu_c, gamma)*chi( dimlesstheta, gamma)* pow(fabs(Z_sec_der(nu, nu_c, gamma)),(-1.0/2.0))
    }
    else
    {
        val=(pow(M_PI,(3.0/2.0))*pow(CHARGE_EL, 2)/(pow(2,(3.0/2.0))*C_LIGHT))*sqrt(nu*nu_c)*n_el_MJ(el_dens, dimlesstheta, gamma)* Z(nu, nu_c, gamma)*chi( dimlesstheta, gamma)* pow(fabs(Z_sec_der(nu, nu_c, gamma)),(-1.0/2.0))
    }
    
    return val;
}

//the functions here are to calculate the total absorption cross section from Ghisellini+ 1991
double C(double nu_ph, double nu_c, double gamma_el, double p_el)
{
    return ((2.0*pow(gamma_el,2)-1)/(gamma_el*pow(p_el,2)))+2*nu_ph*((gamma_el/pow(p_el,2))-gamma_el*log((gamma_el+1)/p_el))/nu_c;
}

double G(double gamma_el, double p_el)
{
    return sqrt(1-2*pow(p_el,2)*(gamma_el*log((gamma_el+1)/p_el)-1));
}

double G_prime(double gamma_el, double p_el)
{
    return (3*gamma_el-(3*pow(gamma_el,2)-1)*log((gamma_el+1)/p_el))/G(gamma_el, p_el)
}

