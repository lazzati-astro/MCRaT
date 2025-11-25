//
// Created by Tyler Parsotan on 11/10/25.
//

#include "mcrat.h"

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr)
{
    //generates an electron with random energy
    double gamma=0, beta=0, phi=0, theta=0, ph_theta=0, ph_phi=0;
    gsl_matrix *rot= gsl_matrix_calloc (3, 3); //create matrix thats 3x3 to do rotation
    gsl_vector_view el_p_prime ; //create vector to hold rotated electron 4 momentum
    gsl_vector *result=gsl_vector_alloc (3);

    gamma=sampleThermalElectron(temp, rand, fPtr);

    //fprintf(fPtr,"Chosen Gamma: %e\n",gamma);

    beta=sqrt( 1- (1/(gamma*gamma)) );
    //printf("Beta is: %e in singleElectron\n", beta);
    phi=gsl_rng_uniform(rand)*2*M_PI;

    theta=sampleElectronTheta(beta, rand, fPtr);
    //fprintf(fPtr,"Beta: %e\tPhi: %e\tTheta: %e\n",beta,phi, theta);
    //fill in electron 4 momentum NOT SURE WHY THE ORDER IS AS SUCH SEEMS TO BE E/c, pz,py,px!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    *(el_p+0)=gamma*(M_EL)*(C_LIGHT);
    *(el_p+1)=gamma*(M_EL)*(C_LIGHT)*beta*cos(theta);
    *(el_p+2)=gamma*(M_EL)*(C_LIGHT)*beta*sin(theta)*sin(phi);
    *(el_p+3)=gamma*(M_EL)*(C_LIGHT)*beta*sin(theta)*cos(phi);

    //printf("Old: %e, %e, %e,%e\n", *(el_p+0), *(el_p+1), *(el_p+2), *(el_p+3));

    el_p_prime=gsl_vector_view_array((el_p+1), 3);

    //find angles of photon NOT SURE WHY WERE CHANGING REFERENCE FRAMES HERE???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ph_phi=atan2(*(ph_p+2), *(ph_p+3)); //Double Check
    ph_theta=atan2(sqrt( pow(*(ph_p+2),2)+  pow(*(ph_p+3),2)) , (*(ph_p+1)) );

    //printf("Calculated Photon phi and theta in singleElectron:%e, %e\n", ph_phi, ph_theta);

    //fill in rotation matrix to rotate around x axis to get rid of phi angle
    gsl_matrix_set(rot, 1,1,1);
    gsl_matrix_set(rot, 2,2,cos(ph_theta));
    gsl_matrix_set(rot, 0,0,cos(ph_theta));
    gsl_matrix_set(rot, 0,2,-sin(ph_theta));
    gsl_matrix_set(rot, 2,0,sin(ph_theta));
    gsl_blas_dgemv(CblasNoTrans, 1, rot, &el_p_prime.vector, 0, result);

    /*
    printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot, 0,0), gsl_matrix_get(rot, 0,1), gsl_matrix_get(rot, 0,2));
    printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot, 1,0), gsl_matrix_get(rot, 1,1), gsl_matrix_get(rot, 1,2));
    printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot, 2,0), gsl_matrix_get(rot, 2,1), gsl_matrix_get(rot, 2,2));

    printf("Middle: %e, %e, %e,%e\n", *(el_p+0), gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2));
    */

    gsl_matrix_set_all(rot,0);

    gsl_matrix_set(rot, 0,0,1);
    gsl_matrix_set(rot, 1,1,cos(-ph_phi));
    gsl_matrix_set(rot, 2,2,cos(-ph_phi));
    gsl_matrix_set(rot, 1,2,-sin(-ph_phi));
    gsl_matrix_set(rot, 2,1,sin(-ph_phi));
    gsl_blas_dgemv(CblasNoTrans, 1, rot, result, 0, &el_p_prime.vector);
    /*
    printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot, 0,0), gsl_matrix_get(rot, 0,1), gsl_matrix_get(rot, 0,2));
    printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot, 1,0), gsl_matrix_get(rot, 1,1), gsl_matrix_get(rot, 1,2));
    printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot, 2,0), gsl_matrix_get(rot, 2,1), gsl_matrix_get(rot, 2,2));
    printf("Final EL_P_vec: %e, %e, %e,%e\n", *(el_p+0), gsl_vector_get(&el_p_prime.vector,0), gsl_vector_get(&el_p_prime.vector,1), gsl_vector_get(&el_p_prime.vector,2));
    */


    gsl_matrix_free (rot);gsl_vector_free(result);
}

double sampleElectronTheta(double beta, gsl_rng * rand, FILE *fPtr)
{
    double y_dum=0, f_x_dum=0, x_dum=0, beta_x_dum=0, theta=0;

    /*this loop is inefficient
    y_dum=1; //initalize loop to get a random theta
    f_x_dum=0;
    while (y_dum>f_x_dum)
    {
        y_dum=gsl_rng_uniform(rand)*1.3;
        x_dum=gsl_rng_uniform(rand)*M_PI;
        f_x_dum=sin(x_dum)*(1-(beta*cos(x_dum)));
    }
    theta=x_dum;
    */

    //can change to this: equation 56 of the RAIKOU paper: DOI: 10.3847/1538-4357/acc94a
    theta = acos((1-sqrt(1+beta*beta+2*beta-4*beta*gsl_rng_uniform(rand)))/beta);


    return theta;
}

double sampleThermalElectron(double temp, gsl_rng * rand, FILE *fPtr)
{
    double gamma=1, factor=0, x_dum=0, y_dum=0, f_x_dum=0, beta_x_dum=0;

    //fprintf(fPtr, "Temp in singleElectron: %e\n", temp);
    if (temp>= 1e7)
    {
        // see also rejection sampling method here: https://arxiv.org/pdf/2408.09105
        //printf("In if\n");
        factor=K_B*temp/(M_EL*C_LIGHT*C_LIGHT);
        y_dum=1; //initalize loop to get a random gamma from the distribution of electron velocities
        f_x_dum=0;
        while ((isnan(f_x_dum) !=0) || (y_dum>f_x_dum) )
        {

            x_dum=gsl_rng_uniform_pos(rand)*(1+100*factor);
            beta_x_dum=sqrt(1-(1/(x_dum*x_dum)));
            y_dum=gsl_rng_uniform(rand)/2.0;

            f_x_dum=x_dum*x_dum*(beta_x_dum/gsl_sf_bessel_Kn (2, 1.0/factor))*exp(-1*x_dum/factor); //
            //fprintf(fPtr,"Choosing a Gamma: xdum: %e, f_x_dum: %e, y_dum: %e\n", x_dum, f_x_dum, y_dum);
        }
        gamma=x_dum;

    }
    else
    {

        //printf("In else\n");
        factor=sqrt(K_B*temp/M_EL);
        //calculate a random gamma from 3 random velocities drawn from a gaussian distribution with std deviation of "factor"
        gamma=1.0/sqrt( 1- (pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+ pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)+pow(gsl_ran_gaussian(rand, factor)/C_LIGHT, 2)  )); //each vel direction is normal distribution -> maxwellian when multiplied
    }

    return gamma;
}

double sampleNonthermalElectron(double p, gsl_rng * rand, FILE *fPtr)
{

    return 0;
}

double samplePowerLaw(double p, double gamma_min, double gamma_max, gsl_rng * rand, FILE *fPtr)
{
    // p: power-law index, gmin/gmax: min/max gamma
    double random_num = gsl_rng_uniform_pos(rand);
    double gamma_e=0;

    if (fabs(p-1.0) < 1e-6)
    {
        // Special case: p ~ 1
        gamma_e = gamma_min * pow(gamma_max / gamma_min, random_num);
    }
    else
    {
        gamma_e = 1.0 + random_num * (pow(gamma_max/gamma_min, 1.0 - p) - 1.0);
        gamma_e = gamma_min * pow(gamma_e, 1.0 / (1.0 - p));
    }
    return gamma_e;
}

double sampleBrokenPowerLaw(double p1, double p2, double gamma_min, double gamma_max, double gamma_break, gsl_rng * rand, FILE *fPtr)
{
    bool p1_is_1=(fabs(p1-1.0) < 1e-6), p2_is_1=(fabs(p2-1.0) < 1e-6);
    double xi_break=0, A=0, gamma_e=0;
    double random_num=gsl_rng_uniform_pos(rand);

    if (!p1_is_1 && !p2_is_1)
    {
        A = 1.0 / ( ((pow(gamma_break, 1.0-p1)-pow(gamma_min, 1.0-p1))/(1.0-p1)) + pow(gamma_break, -1.0*p1+p2)* ((pow(gamma_max, 1.0-p2)-pow(gamma_break, 1.0-p2))/(1.0-p2)));
        xi_break=A*((pow(gamma_break, 1.0-p1)-pow(gamma_min, 1.0-p1))/(1.0-p1));

        if (random_num <= xi_break)
        {
            gamma_e = pow(pow(gamma_min, 1.0-p1) + ((1.0-p1)*random_num/A), 1.0/(1.0-p1));
        }
        else
        {
            // in raikou paper we have: pow(gamma_break, 1.0-p2) - (1.0-p2)...
            // This wasnt sampling the distribution properly and when I set p1=p2=p and tried to recover the single power law
            // behavior found that there should be a plus sign instead to recover the limiting case
            gamma_e = pow( pow(gamma_break, 1.0-p2) + (1.0-p2)*pow(gamma_break, p1-p2)*( ((pow(gamma_min, 1.0-p1)-pow(gamma_break, 1.0-p1))/(1.0-p1))  +(random_num/A))  , 1.0/(1.0-p2));
        }
    }
    else if (p1_is_1 && !p2_is_1)
    {
        A = 1.0 / (log(gamma_break/gamma_min) + pow(gamma_break, -1.0*p1+p2)*((pow(gamma_max, 1.0-p2)-pow(gamma_break, 1.0-p2))/(1.0-p2)) ) ;
        xi_break=A*log(gamma_break/gamma_min);

        if (random_num <= xi_break)
        {
            gamma_e = gamma_min * exp(random_num/A);
        }
        else
        {
            gamma_e = pow( pow(gamma_break, 1.0-p2) - (1.0-p2)*pow(gamma_break, p1-p2)*(log(gamma_break/gamma_min) -(random_num/A))  , 1.0/(1.0-p2));
        }

    }
    else if (!p1_is_1 && p2_is_1)
    {
        A = 1.0/(((pow(gamma_break, 1.0-p1)-pow(gamma_min, 1.0-p1))/(1.0-p1)) + pow(gamma_break, -1.0*p1+p2)*log(gamma_max/gamma_break));
        xi_break = A * ((pow(gamma_break, 1.0-p1)-pow(gamma_min, 1.0-p1))/(1.0-p1));

        if (random_num <= xi_break)
        {
            gamma_e = pow(pow(gamma_min, 1.0-p1) + ((1.0-p1)*random_num/A), 1.0/(1.0-p1));
        }
        else
        {
            gamma_e = gamma_break*exp(pow(gamma_break, p1-p2)* ((random_num/A) - ((pow(gamma_break, 1.0-p1)-pow(gamma_min, 1.0-p1))/(1.0-p1)) ) );
        }
    }
    else
    {
        fprintf(fPtr,"sampleBrokenPowerLaw: In the else\n p1_is_1: %d, p2_is_1: %d \n This shouldnt have occured exiting", p1_is_1, p2_is_1);
        exit(1);
    }

    return gamma_e;

}

double brokenPowerLawNorm(double p1, double p2, double gamma_min, double gamma_max, double gamma_break)
{
    bool p1_is_1 = (fabs(p1 - 1.0) < 1e-10);
    bool p2_is_1 = (fabs(p2 - 1.0) < 1e-10);
    double A;

    // Case 1: Neither p1 nor p2 equals 1
    if (!p1_is_1 && !p2_is_1)
    {
        double term1 = (pow(gamma_break, 1.0 - p1) - pow(gamma_min, 1.0 - p1)) / (1.0 - p1);
        double term2 = pow(gamma_break, -p1 + p2) *
                       (pow(gamma_max, 1.0 - p2) - pow(gamma_break, 1.0 - p2)) / (1.0 - p2);
        A = 1.0 / (term1 + term2);
    }
    // Case 2: p1 equals 1, p2 does not
    else if (p1_is_1 && !p2_is_1)
    {
        double term1 = log(gamma_break / gamma_min);
        double term2 = pow(gamma_break, p2 - 1.0) *
                       (pow(gamma_max, 1.0 - p2) - pow(gamma_break, 1.0 - p2)) / (1.0 - p2);
        A = 1.0 / (term1 + term2);
    }
    // Case 3: p1 does not equal 1, p2 equals 1
    else if (!p1_is_1 && p2_is_1)
    {
        double term1 = (pow(gamma_break, -p1 + 1.0) - pow(gamma_min, -p1 + 1.0)) / (-p1 + 1.0);
        double term2 = pow(gamma_break, -p1 + 1.0) * log(gamma_max / gamma_break);
        A = 1.0 / (term1 + term2);
    }
    // Case 4: Both equal 1 (should not occur in practice)
    else
    {
        // Handle this edge case - could use double logarithmic form
        A = 0;
    }

    return A;
}


double singleElectronBrokenPowerLaw(double x, double p1, double p2, double gamma_min, double gamma_max, double gamma_break)
{
    /**
     * Calculate broken power-law distribution value for a single x value
     *
     * @param x Gamma value (Lorentz factor)
     * @param p1 Power-law index below break
     * @param p2 Power-law index above break
     * @param gamma_min Minimum Lorentz factor
     * @param gamma_max Maximum Lorentz factor
     * @param gamma_break Break Lorentz factor
     * @return Distribution value f(x)
     */

    // Calculate normalization constant
    double A = brokenPowerLawNorm(p1, p2, gamma_min, gamma_max, gamma_break);

    // Calculate distribution value based on which segment x falls in
    double y;
    if (x <= gamma_break)
    {
        // Low-energy segment
        y = A * pow(x, -p1);
    }
    else
    {
        // High-energy segment with continuity factor
        double continuity_factor = pow(gamma_break, -p1 + p2);
        y = A * pow(x, -p2) * continuity_factor;
    }

    return y;
}


void arrayElectronBrokenPowerLaw(const double *x, double *y, int n_points, double p1, double p2, double gamma_min, double gamma_max, double gamma_break)
{
    /**
     * Calculate broken power-law distribution values for an array of x values
     *
     * @param x Array of gamma values (input, must be pre-allocated)
     * @param y Array of distribution values (output, must be pre-allocated)
     * @param n_points Number of points in x and y arrays
     * @param p1 Power-law index below break
     * @param p2 Power-law index above break
     * @param gamma_min Minimum Lorentz factor
     * @param gamma_max Maximum Lorentz factor
     * @param gamma_break Break Lorentz factor
     */

    // Calculate normalization constant once (more efficient)
    double A = brokenPowerLawNorm(p1, p2, gamma_min, gamma_max, gamma_break);

    // Calculate continuity factor once
    double continuity_factor = pow(gamma_break, -p1 + p2);

    // Fill the distribution values
    for (int i = 0; i < n_points; i++)
    {
        if (x[i] <= gamma_break)
        {
            // Low-energy segment
            y[i] = A * pow(x[i], -p1);
        }
        else
        {
            // High-energy segment
            y[i] = A * pow(x[i], -p2) * continuity_factor;
        }
    }
}


double powerLawNorm(double p, double gamma_min, double gamma_max)
{
    /**
     * Calculate the normalization constant A for a single power-law distribution
     *
     * The distribution is: n(gamma) = A * gamma^(-p) for gamma_min <= gamma <= gamma_max
     *
     * @param p Power-law index
     * @param gamma_min Minimum Lorentz factor
     * @param gamma_max Maximum Lorentz factor
     * @return Normalization constant A
     */

    bool p_is_1 = (fabs(p - 1.0) < 1e-10);
    double A;

    if (p_is_1)
    {
        // Special case: p = 1 (logarithmic distribution)
        A = 1.0 / log(gamma_max / gamma_min);
    }
    else
    {
        // General case: p != 1
        double integral = (pow(gamma_max, 1.0 - p) - pow(gamma_min, 1.0 - p)) / (1.0 - p);
        A = 1.0 / integral;
    }

    return A;
}


double singleElectronPowerLaw(double x, double p, double gamma_min, double gamma_max)
{
    /**
     * Calculate single power-law distribution value for a single x value
     *
     * @param x Gamma value (Lorentz factor)
     * @param p Power-law index
     * @param gamma_min Minimum Lorentz factor
     * @param gamma_max Maximum Lorentz factor
     * @return Distribution value f(x)
     */

    // Check if x is within valid range
    if (x < gamma_min || x > gamma_max)
    {
        return 0.0;  // Outside valid range
    }

    // Calculate normalization constant
    double A = powerLawNorm(p, gamma_min, gamma_max);

    // Calculate distribution value
    double y = A * pow(x, -p);

    return y;
}



void arrayElectronPowerLaw(const double *x, double *y, int n_points, double p, double gamma_min, double gamma_max)
{
    /**
     * Calculate single power-law distribution values for an array of x values
     *
     * @param x Array of gamma values (input, must be pre-allocated)
     * @param y Array of distribution values (output, must be pre-allocated)
     * @param n_points Number of points in x and y arrays
     * @param p Power-law index
     * @param gamma_min Minimum Lorentz factor
     * @param gamma_max Maximum Lorentz factor
     */

    // Calculate normalization constant once (more efficient)
    double A = powerLawNorm(p, gamma_min, gamma_max);

    // Fill the distribution values
    for (int i = 0; i < n_points; i++)
    {
        if (x[i] >= gamma_min && x[i] <= gamma_max)
        {
            y[i] = A * pow(x[i], -p);
        }
        else
           {
            y[i] = 0.0;  // Outside valid range
        }
    }
}

double singleMaxwellJuttner(double gamma, double theta)
{
    /*
       returns the maxwell juttner distribution value for a given dimensionless temp, theta, and electron lorentz factor, gamma
        When the temp >~6e7 K, we use the full distribution since we evaluate the integral at this threshold and it is
        properly normalized (there is a ~12% error in the normalization at this threshold value). Below this temperature, we calculated the
        limit of the distribution as theta->0 and the normalizations
    */
    double result=0, normalization=0;
    
    if (theta > 1.e-2) 
    {
        normalization = gsl_sf_bessel_Kn(2, 1. / theta) * exp(1. / theta);
    } 
    else 
    {
        normalization = sqrt(M_PI * theta / 2.);
    }

    result = ((gamma * sqrt(gamma * gamma - 1.) / (theta * normalization)) * exp(-(gamma - 1.) / theta));

    return result;
}

double nonThermalElectronDistIntegrand(double x, void * params)
{
    double result=0;
    #if NONTHERMAL_E_DIST != OFF

        #if NONTHERMAL_E_DIST == POWERLAW
            result = singleElectronPowerLaw(x, POWERLAW_INDEX, GAMMA_MIN, GAMMA_MAX);
        #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
            result = singleElectronBrokenPowerLaw(x, POWERLAW_INDEX_1, POWERLAW_INDEX_2, GAMMA_MIN, GAMMA_MAX, GAMMA_BREAK);
        #else
            #error Unknown nonthermal electron distribution.
        #endif
    #endif

    return result;
}



double calculateNormPowerLawEnergyDens(double p, double gamma_min, double gamma_max)
{
    //calcualte \gamma*m*c^2 times the power law electron distribution
    double result=0;
    bool p_is_2 = (fabs(p - 2.0) < 1e-10);

    if (p_is_2)
    {
        // Special case: p = 2 (logarithmic distribution)
        result = log(gamma_max / gamma_min);
    }
    else
    {
        // General case: p != 2
        result = (pow(gamma_max, 2.0 - p) - pow(gamma_min, 2.0 - p)) / (2.0 - p);
    }

    //now multiply by the actual normalization of the distribution
    result*=powerLawNorm(p, gamma_min, gamma_max);

    //now multiply by electron rest mass
    result*=(M_EL*C_LIGHT*C_LIGHT);
    return  result;

}

double calculateNormBrokenPowerLawEnergyDens(double p1, double p2, double gamma_min, double gamma_max, double gamma_break)
{
    //calcualte \gamma*m*c^2 times the broken power law electron distribution
    bool p1_is_2 = (fabs(p1 - 2.0) < 1e-10);
    bool p2_is_2 = (fabs(p2 - 2.0) < 1e-10);
    double result=0, term1, term2;

    // Case 1: Neither p1 nor p2 equals 1
    if (!p1_is_2 && !p2_is_2)
    {
        term1 = (pow(gamma_break, 2.0 - p1) - pow(gamma_min, 2.0 - p1)) / (2.0 - p1);
        term2 = pow(gamma_break, -p1 + p2) *
                       (pow(gamma_max, 2.0 - p2) - pow(gamma_break, 2.0 - p2)) / (2.0 - p2);
        result=term1 + term2;
    }
    // Case 2: p1 equals 1, p2 does not
    else if (p1_is_2 && !p2_is_2)
    {
        term1 = log(gamma_break / gamma_min);
        term2 = pow(gamma_break, p2 - 2.0) *
                       (pow(gamma_max, 2.0 - p2) - pow(gamma_break, 2.0 - p2)) / (2.0 - p2);
        result = term1 + term2;
    }
    // Case 3: p1 does not equal 1, p2 equals 1
    else if (!p1_is_2 && p2_is_2)
    {
        term1 = (pow(gamma_break, -p1 + 2.0) - pow(gamma_min, -p1 + 2.0)) / (-p1 + 2.0);
        term2 = pow(gamma_break, -p1 + 2.0) * log(gamma_max / gamma_break);
        result = (term1 + term2);
    }
    // Case 4: Both equal 2 (should not occur in practice)
    else
    {
        // Handle this edge case - could use double logarithmic form
        result = 0;
    }

    //now multiply by the actual normalization of the distribution
    result*=brokenPowerLawNorm(p1, p2, gamma_min, gamma_max, gamma_break);

    //now multiply by electron rest mass
    result*=(M_EL*C_LIGHT*C_LIGHT);

    return result;

}

#if NONTHERMAL_E_DIST != OFF
    void calculateElectronDistSubgroupDens(double *subgroup_dens, FILE *fPtr)
    {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        int i;
        double result, error;
        double dgamma = (log10(GAMMA_MAX) - log10(GAMMA_MIN)) / N_GAMMA;
        double gamma_min, gamma_max;
        gsl_function F;
        F.function = &nonThermalElectronDistIntegrand;

        for (i=0;i<N_GAMMA; i++)
        {
            gamma_min = pow(10.0, log10(GAMMA_MIN) + i * dgamma);
            gamma_max = pow(10.0, log10(GAMMA_MIN) + (i+1) * dgamma);

            gsl_integration_qags (&F, gamma_min, gamma_max, 0, 1e-7, 1000, w, &result, &error);
            fprintf(fPtr, "gamma_min: %e gamma_max: %e, result: %e\n", gamma_min, gamma_max, result);

            subgroup_dens[i] = result;
        }
    }

    void calculateNonthermalElectronDens(struct hydro_dataframe *hydro_data, FILE *fPtr)
    {
        /*
            set the non-thermal electron energy density to the magnetic field energy density
        */
        int i=0;
        double result=0;
        double b_field = 0;
        double energy_dens_per_particle=0;

        for (i=0; i<hydro_data->num_elements; i++)
        {

            b_field = getMagneticFieldMagnitude(hydro_data, i)
            #if NONTHERMAL_E_DIST == POWERLAW
                energy_dens_per_particle = calculateNormPowerLawEnergyDens(POWERLAW_INDEX, GAMMA_MIN, GAMMA_MAX);
            #elif NONTHERMAL_E_DIST == BROKENPOWERLAW
                energy_dens_per_particle = calculateNormBrokenPowerLawEnergyDens(POWERLAW_INDEX_1, POWERLAW_INDEX_2, GAMMA_MIN, GAMMA_MAX, GAMMA_BREAK);
            #else
                #error Unknown nonthermal electron distribution.
            #endif

            //this is the number density
            (hydro_data->nonthermal_dens)[i] =  (b_field * b_field)/(8.0*M_PI* energy_dens_per_particle);

            fprintf("in cell %d the particle number density is: %e, the nonthermal number density is: %e\n", i, (hydro_data->dens)[i]/M_P, (hydro_data->nonthermal_dens)[i]);
        }
    }

#endif
