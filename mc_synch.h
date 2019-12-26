/*
This header file is for the different functions for emitting and absorbing synchrotron photons
*/

double calcCyclotronFreq(double magnetic_field);

double calcDimlessTheta(double temp);

double calcB(double el_dens, double temp, double epsilon_b);

double n_el_MJ(double el_dens, double dimlesstheta, double gamma);

double n_el_MB(double el_dens, double dimlesstheta, double gamma);

double Z(double nu, double nu_c, double gamma );

double Z_sec_der(double nu, double nu_c, double gamma);

double chi(double dimlesstheta, double gamma);

double gamma0(double nu, double nu_c, double dimlesstheta);

double jnu(double nu, double nu_c, double dimlesstheta, double el_dens);
