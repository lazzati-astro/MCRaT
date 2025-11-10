//
// Created by Tyler Parsotan on 11/10/25.
//

#include "mcrat.h"

void singleElectron(double *el_p, double temp, double *ph_p, gsl_rng * rand, FILE *fPtr)
{
    //generates an electron with random energy
    double factor=0, gamma=0;
    double y_dum=0, f_x_dum=0, x_dum=0, beta_x_dum=0, beta=0, phi=0, theta=0, ph_theta=0, ph_phi=0;
    gsl_matrix *rot= gsl_matrix_calloc (3, 3); //create matrix thats 3x3 to do rotation
    gsl_vector_view el_p_prime ; //create vector to hold rotated electron 4 momentum
    gsl_vector *result=gsl_vector_alloc (3);

    gamma=sampleThermalElectron(temp, rand, fPtr);

    //fprintf(fPtr,"Chosen Gamma: %e\n",gamma);

    beta=sqrt( 1- (1/(gamma*gamma)) );
    //printf("Beta is: %e in singleElectron\n", beta);
    phi=gsl_rng_uniform(rand)*2*M_PI;

    y_dum=1; //initalize loop to get a random theta
    f_x_dum=0;
    while (y_dum>f_x_dum)
    {
        y_dum=gsl_rng_uniform(rand)*1.3;
        x_dum=gsl_rng_uniform(rand)*M_PI;
        f_x_dum=sin(x_dum)*(1-(beta*cos(x_dum)));
    }
    theta=x_dum;
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