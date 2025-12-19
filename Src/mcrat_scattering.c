//
//  mcrat_scattering.c
//  file contains all functions for MCRaT to conduct a compton/Klein-nishina scattering including polarization
//
//  Created by Tyler Parsotan on 7/23/21.
//

#include "mcrat.h"

void mullerMatrixRotation(double theta, double *s, FILE *fPtr)
{
    //makes a CCW rotation od the stokes parameters when the photon velocity vector is pointed towards the observer, follows Lundman
    gsl_matrix *M= gsl_matrix_calloc (4, 4); //create matrix thats 4x4 to do rotation as defined in McMaster 1961 (has it to rotate CW in that paper)
    gsl_vector *result= gsl_vector_calloc(4);
    gsl_vector_view stokes;
    
    stokes=gsl_vector_view_array(s, 4);
    //fprintf(fPtr, "sokes parameter before= %e %e %e %e\n", gsl_vector_get(&stokes.vector, 0), gsl_vector_get(&stokes.vector, 1), gsl_vector_get(&stokes.vector, 2), gsl_vector_get(&stokes.vector, 3));
    
    gsl_matrix_set(M, 0,0,1);
    gsl_matrix_set(M, 3,3,1);
    gsl_matrix_set(M, 1,1,cos(2*theta));
    gsl_matrix_set(M, 2,2,cos(2*theta));
    gsl_matrix_set(M, 1,2,-1*sin(2*theta));
    gsl_matrix_set(M, 2,1,sin(2*theta));
    gsl_blas_dgemv(CblasNoTrans, 1, M, &stokes.vector, 0, result); //Ms=s
    
    //fprintf(fPtr, "stokes parameter after= %e %e %e %e\n\n", gsl_vector_get(result, 0), gsl_vector_get(result, 1), gsl_vector_get(result, 2), gsl_vector_get(result, 3));
    
    //save back to the original stokes vector
    *(s+0)=gsl_vector_get(result, 0);
    *(s+1)=gsl_vector_get(result, 1);
    *(s+2)=gsl_vector_get(result, 2);
    *(s+3)=gsl_vector_get(result, 3);
    
    gsl_vector_free(result);
    gsl_matrix_free (M);
    
}

void findXY(double *v_ph, double *vector, double *x, double *y)
{
    //finds the stokes plane coordinate x,y axis for the photon velocity with respect to some reference vector
    //assumes that pointers point to array of 3 doubles in length
    double norm=0;
    
    *(y+0)= ((*(v_ph+1))*(*(vector+2))-(*(v_ph+2))*(*(vector+1)));
    *(y+1)= -1*((*(v_ph+0))*(*(vector+2))-(*(v_ph+2))*(*(vector+0)));
    *(y+2)= ((*(v_ph+0))*(*(vector+1))-(*(v_ph+1))*(*(vector+0))); // vector X v_ph
    
    if ( (*(y+0)==0) && (*(y+1)==0)|| (*(y+2))==0)
    {
        printf("The calculated y stokes plane coordinate value is 0\n\n");
        printf("This is most likely due to the boosted photon velocity vector being parallel to the velocity vector that was used for the lorentz boost. \n\n");
    }

    
    norm=1.0/sqrt( (*(y+0))*(*(y+0)) + (*(y+1))*(*(y+1)) + (*(y+2))*(*(y+2)));
    *(y+0) *= norm;
    *(y+1) *= norm;
    *(y+2) *= norm;
    
    *(x+0)= (*(y+1))*(*(v_ph+2))-(*(y+2))*(*(v_ph+1));
    *(x+1)= -1*((*(y+0))*(*(v_ph+2))-(*(y+2))*(*(v_ph+0)));
    *(x+2)= (*(y+0))*(*(v_ph+1))-(*(y+1))*(*(v_ph+0));
    
    norm=1.0/sqrt( (*(x+0))*(*(x+0)) + (*(x+1))*(*(x+1)) + (*(x+2))*(*(x+2)));
    *(x+0) *= norm;
    *(x+1) *= norm;
    *(x+2) *= norm;
    
}

double findPhi(double *x_old, double *y_old, double *x_new, double *y_new)
{
    //find the angle to rotate the stokes vector to transform from one set of stokes coordinates to another
    //this is given by Lundman
    gsl_vector_view y=gsl_vector_view_array(y_old, 3);
    gsl_vector_view x=gsl_vector_view_array(x_old, 3);
    gsl_vector_view y_prime=gsl_vector_view_array(y_new, 3);
    gsl_vector_view x_prime=gsl_vector_view_array(x_new, 3);
    double factor=0, dot_prod_result=0;
    
    gsl_blas_ddot(&x.vector, &y_prime.vector, &dot_prod_result);
    
    if (dot_prod_result>0)
    {
        factor=1;
    }
    else if (dot_prod_result<0)
    {
        factor=-1;
    }
    else
    {
        factor=0;
    }
    
    gsl_blas_ddot(&y.vector, &y_prime.vector, &dot_prod_result);
    
    if ((dot_prod_result<-1) || (dot_prod_result>1))
    {
        //printf("The old dot poduct was %e, the new one is %e\n",dot_prod_result, round(dot_prod_result));
        dot_prod_result=round(dot_prod_result);//do this rounding so numerical error that causes value to be <-1 or >1 gets rounded and becomes a real value if its close enough to these limits
    }
    
    return -1*factor*acos(dot_prod_result);
}

void stokesRotation(double *v, double *v_ph, double *v_ph_boosted, double *s, FILE *fPtr)
{
    //takes 3 velocities of the initial photon, v_ph, the boosted photon, v_ph_boosted. and the boost vector, v
    double z_hat[3]={0,0,1}; //z to calulate stokes
    double x[3]={0,0,0}, y[3]={0,0,0}, x_new[3]={0,0,0}, y_new[3]={0,0,0};//initalize arrays to hold stokes coordinate system
    double phi=0;
    
    //if (i==0)
    {
    //find stokes coordinate sys in orig frame with respect to z axis
    findXY(v_ph, &z_hat, &x, &y);
    }
    
    //find stokes coordinate sys in orig frame with respect to boost vector
    findXY(v_ph, v, &x_new, &y_new);
    
    phi=findPhi(x, y, x_new, y_new);//now find rotation between the two coordinate systems
    
    //rotate the stokes vector now to put it in the coordinate system fo the boosted photon and the boost evctor
    mullerMatrixRotation(phi, s, fPtr);
    
    
    if ( isnan(*(s+0)) || isnan(*(s+1)) || isnan(*(s+2)) || isnan(*(s+3)) )
    {
        printf("A stokes value is nan\n\n");
    }
     
    
    //find the new coordinates of the rotated stokes vector with the boosted photon and the boost vector
    findXY(v_ph_boosted, v, &x, &y);
    
    if ( isnan(*(y+0)) || isnan(*(y+1)) || isnan(*(y+2)))
    {
        printf("A plane coordinate value is nan\n\n");
        printf("This is most likely due to the boosted photon velocity vector being parallel to the velocity vector that was used for the lorentz boost. \n\n");
        fprintf(fPtr, "A plane coordinate value is nan\n\n");
        fprintf(fPtr, "This is most likely due to the boosted photon velocity vector being parallel to the velocity vector that was used for the lorentz boost. \n\n");
        fflush(fPtr);
        memcpy(y, y_new, 3*sizeof(double));
        memcpy(x, x_new, 3*sizeof(double));
    }

    
    //find stokes coordinate sys in orig frame with respect to z axis
    findXY(v_ph_boosted, &z_hat, &x_new, &y_new);
    
    phi=findPhi(x, y, x_new, y_new);//now find rotation between the two coordinate systems
    
    //do the rotation of the stokes vector to put it in the coordinate system of the boosted photon and the z axis
    mullerMatrixRotation(phi, s, fPtr);
    
    /*
    if ( isnan(*(s+0)) || isnan(*(s+1)) || isnan(*(s+2)) || isnan(*(s+3)) )
    {
        printf("A stokes value is nan\n\n");
    }
     */
    
}

void stokesScatter(double *s, double *orig_s,  gsl_vector *ph_p_orig, gsl_vector *result0, double ph_p_prime, double scattered_ph_e, FILE *fPtr)
{
    double dot_prod_result=0;
    double *z_axis_electron_rest_frame=malloc(3*sizeof(double)); //basis vector of the z axis in the elctron rest frame
    double x_tilde[3]={0,0,0}, y_tilde[3]={0,0,0}, x_tilde_new[3]={0,0,0}, y_tilde_new[3]={0,0,0};//initalize arrays to hold stokes coordinate system
    gsl_matrix *scatt= gsl_matrix_calloc (4, 4); //fano's matrix for scattering stokes parameters
    gsl_vector *scatt_result=gsl_vector_calloc (4);
    gsl_vector_view stokes;


    stokes=gsl_vector_view_array(orig_s, 4);

    
    //fill in z-axis basis vector
    *(z_axis_electron_rest_frame+0)=0;
    *(z_axis_electron_rest_frame+1)=0;
    *(z_axis_electron_rest_frame+2)=1;

    //orient the stokes coordinate system such that its perpendicular to the scattering plane
    findXY(gsl_vector_ptr(ph_p_orig, 1),z_axis_electron_rest_frame, x_tilde, y_tilde);
    findXY(gsl_vector_ptr(result0,0),gsl_vector_ptr(ph_p_orig, 1), x_tilde_new, y_tilde_new);
    phi=findPhi(x_tilde, y_tilde, x_tilde_new, y_tilde_new);
    mullerMatrixRotation(phi, s, fPtr);
    
    //find the theta between the incoming and scattered photons, by doing dot product and taking arccos of it
    dot_prod_result=(gsl_vector_get(ph_p_orig,1)*gsl_vector_get(result0,0)+gsl_vector_get(ph_p_orig,2)*gsl_vector_get(result0,1)+gsl_vector_get(ph_p_orig,3)*gsl_vector_get(result0,2) )/(gsl_vector_get(ph_p_orig,0)*(ph_p_prime)) ;
    
    if ((dot_prod_result<-1) || (dot_prod_result>1))
    {
        //printf("The old dot poduct was %e, the new one is %e\n",dot_prod_result, round(dot_prod_result));
        dot_prod_result=round(dot_prod_result);//do this rounding so numerical error that causes value to be <-1 or >1 gets rounded and becomes a real value if its close enough to these limits
    }
    theta=acos(dot_prod_result);

    
    //do the scattering of the stokes parameters
    gsl_matrix_set(scatt, 0,0,1.0+pow(cos(theta), 2.0)+((1-cos(theta))*(gsl_vector_get(ph_p_orig,0) - gsl_vector_get(result,0))/(M_EL*C_LIGHT ) ) ); //following lundman's matrix
    gsl_matrix_set(scatt, 0,1, sin(theta)*sin(theta));
    gsl_matrix_set(scatt, 1,0, sin(theta)*sin(theta));
    gsl_matrix_set(scatt, 1,1,1.0+cos(theta)*cos(theta));
    gsl_matrix_set(scatt, 2,2, 2.0*cos(theta));
    gsl_matrix_set(scatt, 3,3, 2.0*cos(theta)+ ((cos(theta))*(1-cos(theta))*(gsl_vector_get(ph_p_orig,0) - scattered_ph_e)/(M_EL*C_LIGHT )) );
    //gsl_matrix_scale(scatt, (gsl_vector_get(result,0)/(*(ph_p_prime+0)))*((gsl_vector_get(result,0)/(*(ph_p_prime+0))))*0.5*3*THOM_X_SECT/(8*M_PI) ); //scale the matrix by 0.5*r_0^2 (\epsilon/\epsilon_0)^2 DONT NEED THIS BECAUSE WE NORMALIZE STOKES VECTOR SO THIS CANCELS ITSELF OUT
    gsl_blas_dgemv(CblasNoTrans, 1, scatt, &stokes.vector, 0, scatt_result);
    /*
     fprintf(fPtr,"before s: %e, %e, %e,%e\n", gsl_vector_get(&stokes.vector,0), gsl_vector_get(&stokes.vector,1), gsl_vector_get(&stokes.vector,2), gsl_vector_get(&stokes.vector,3));
     fprintf(fPtr,"Scatt Matrix 0: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 0,0), gsl_matrix_get(scatt, 0,1), gsl_matrix_get(scatt, 0,2), gsl_matrix_get(scatt, 0,3));
     fprintf(fPtr,"Scatt Matrix 1: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 1,0), gsl_matrix_get(scatt, 1,1), gsl_matrix_get(scatt, 1,2), gsl_matrix_get(scatt, 1,3));
     fprintf(fPtr,"Scatt Matrix 2: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 2,0), gsl_matrix_get(scatt, 2,1), gsl_matrix_get(scatt, 2,2), gsl_matrix_get(scatt, 2,3));
     fprintf(fPtr,"Scatt Matrix 3: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 3,0), gsl_matrix_get(scatt, 3,1), gsl_matrix_get(scatt, 3,2), gsl_matrix_get(scatt, 3,3));
     fprintf(fPtr,"s: %e, %e, %e,%e\n", gsl_vector_get(scatt_result,0), gsl_vector_get(scatt_result,1), gsl_vector_get(scatt_result,2), gsl_vector_get(scatt_result,3));
     */
    
    
    //normalize and rotate back
    *(s+0)=gsl_vector_get(scatt_result,0)/gsl_vector_get(scatt_result,0); //should be 1.0
    *(s+1)=gsl_vector_get(scatt_result,1)/gsl_vector_get(scatt_result,0);
    *(s+2)=gsl_vector_get(scatt_result,2)/gsl_vector_get(scatt_result,0);
    *(s+3)=gsl_vector_get(scatt_result,3)/gsl_vector_get(scatt_result,0);
    //fprintf(fPtr,"s after norm: %e, %e, %e,%e\n", gsl_vector_get(&stokes.vector,0), gsl_vector_get(&stokes.vector,1), gsl_vector_get(&stokes.vector,2), gsl_vector_get(&stokes.vector,3));
    
    
    //need to find current stokes coordinate system defined in the plane of k-k_0
    findXY(gsl_vector_ptr(result0,0),gsl_vector_ptr(ph_p_orig, 1), x_tilde, y_tilde);
    
    //then find the new coordinate system between scattered photon 4 onetum and the z axis
    findXY(gsl_vector_ptr(result0,0),z_axis_electron_rest_frame, x_tilde_new, y_tilde_new);
    
    //find phi to transform between the two coodinate systems
    phi=findPhi(x_tilde, y_tilde, x_tilde_new, y_tilde_new);
    
    //do the rotation
    mullerMatrixRotation(phi, s, fPtr);
    
    free(z_axis_electron_rest_frame);
    gsl_matrix_free(scatt);
    gsl_vector_free(scatt_result);
}


int singleScatter(double *el_comov, double *ph_comov, double *s, gsl_rng * rand, FILE *fPtr)
{
    //This routine performs a scattering between a photon and a moving electron.
    int scattering_occured=0;
    double dotprod_1; //to test orthogonality
    double *z_axis_electron_rest_frame=malloc(3*sizeof(double)); //basis vector of the z axis in the elctron rest frame
    double *el_v=malloc(3*sizeof(double));
    double *negative_el_v=malloc(3*sizeof(double));
    double *ph_p_prime=malloc(4*sizeof(double));//use this to keep track of how the ph 4 momentum changes with each rotation
    double *el_p_prime=malloc(4*sizeof(double));
    double phi0=0, phi1=0, phi=0, theta=0;
    double x_tilde[3]={0,0,0}, y_tilde[3]={0,0,0}, x_tilde_new[3]={0,0,0}, y_tilde_new[3]={0,0,0};//initalize arrays to hold stokes coordinate system
    gsl_matrix *rot0= gsl_matrix_calloc (3, 3); //create matricies thats 3x3 to do rotations
    gsl_matrix *rot1= gsl_matrix_calloc (3, 3);
    gsl_matrix *scatt= gsl_matrix_calloc (4, 4); //fano's matrix for scattering stokes parameters
    gsl_vector *scatt_result=gsl_vector_calloc (4);
    gsl_vector *result0=gsl_vector_calloc (3); //vectors to hold results of rotations
    gsl_vector *result1=gsl_vector_calloc (3);
    gsl_vector *result=gsl_vector_calloc (4);
    gsl_vector *whole_ph_p=gsl_vector_calloc (4);
    gsl_vector *ph_p_orig=gsl_vector_calloc (4) ;//vector to hold the original incoming photon velocity vector in the electron rest frame
    gsl_vector_view ph_p ;//create vector to hold comoving photon and electron 4 momentum
    gsl_vector_view el_p ;
    gsl_vector_view stokes, test;
    

    /*
     Dont need these vectors anymore, plus didnt have code to free allocations so it was causing memory leaks
    gsl_vector *result0_x=gsl_vector_alloc (3); //vectors to hold results of rotations for stokes coordinates
    gsl_vector *result1_x=gsl_vector_alloc (3);
    gsl_vector *result0_y=gsl_vector_alloc (3); //vectors to hold results of rotations for stokes coordinates
    gsl_vector *result1_y=gsl_vector_alloc (3);
     */
    
    //fill in z-axis basis vector
    *(z_axis_electron_rest_frame+0)=0;
    *(z_axis_electron_rest_frame+1)=0;
    *(z_axis_electron_rest_frame+2)=1;
    
    /* was for testing against Kraw
    *(s+0)=1; //should be 1.0
    *(s+1)=1;
    *(s+2)=0;
    *(s+3)=0;
    
    *(ph_comov+0)=PL_CONST*1e12/C_LIGHT;
    *(ph_comov+1)=0; //set values of photon prime momentum from doing the scattering to use the vector view of it in dot product
    *(ph_comov+2)=0;
    *(ph_comov+3)=PL_CONST*1e12/C_LIGHT;
    
    theta=85*M_PI/180;
    phi=0;
    dotprod_1=pow(1-(pow(100, -2.0)) ,0.5);
    *(el_comov+0)=100*M_EL*C_LIGHT;
    *(el_comov+1)=100*M_EL*C_LIGHT*dotprod_1*sin(theta)*cos(phi); //set values of photon prime momentum from doing the scattering to use the vector view of it in dot product
    *(el_comov+2)=100*M_EL*C_LIGHT*dotprod_1*sin(theta)*sin(phi);
    *(el_comov+3)=100*M_EL*C_LIGHT*dotprod_1*cos(theta);
     */
    
    //fill in electron velocity array and photon 4 momentum
    *(el_v+0)=(*(el_comov+1))/(*(el_comov+0));
    *(el_v+1)=(*(el_comov+2))/(*(el_comov+0));
    *(el_v+2)=(*(el_comov+3))/(*(el_comov+0));
    //printf("el_v: %e, %e, %e\n", *(el_v+0), *(el_v+1), *(el_v+2));
    
    //lorentz boost into frame where the electron is stationary
    lorentzBoost(el_v, el_comov, el_p_prime, 'e', fPtr);
    lorentzBoost(el_v, ph_comov, ph_p_prime, 'p', fPtr);
    //printf("New ph_p in electron rest frame: %e, %e, %e,%e\n", *(ph_p_prime+0), *(ph_p_prime+1), *(ph_p_prime+2), *(ph_p_prime+3));
    
    //rotate 'stokes plane'
    //if (STOKES_SWITCH != 0)
    #if STOKES_SWITCH == ON
    {
        stokesRotation(el_v, (ph_comov+1), (ph_p_prime+1), s, fPtr);
        stokes=gsl_vector_view_array(s, 4);

    }
    #endif
    
    //printf(fPtr, "y_tilde: %e, %e, %e\n", *(y_tilde+0), *(y_tilde+1), *(y_tilde+2));
    
    ph_p=gsl_vector_view_array((ph_p_prime+1), 3);
    el_p=gsl_vector_view_array(el_p_prime,4);
    
    gsl_vector_set(ph_p_orig, 0, *(ph_p_prime+0));
    gsl_vector_set(ph_p_orig, 1, *(ph_p_prime+1));
    gsl_vector_set(ph_p_orig, 2, *(ph_p_prime+2));
    gsl_vector_set(ph_p_orig, 3, *(ph_p_prime+3));
    
    //gsl_blas_ddot(&y_tilde_rot.vector, &ph_p.vector, &dotprod_1);
    //fprintf(fPtr, "After lorentz boost Angle between the  y_tilde_rot and the photon velocity vector is: %e\n", acos(dotprod_1/ gsl_blas_dnrm2(&ph_p.vector))*180/M_PI);
    
    phi0=atan2(*(ph_p_prime+2), *(ph_p_prime+1) );
    //fprintf(fPtr,"Photon Phi: %e\n", phi0);
    //rotate the axes so that the photon incomes along the x-axis
    gsl_matrix_set(rot0, 2,2,1);
    gsl_matrix_set(rot0, 0,0,cos(-phi0));
    gsl_matrix_set(rot0, 1,1,cos(-phi0));
    gsl_matrix_set(rot0, 0,1,-sin(-phi0));
    gsl_matrix_set(rot0, 1,0,sin(-phi0));
    gsl_blas_dgemv(CblasNoTrans, 1, rot0, &ph_p.vector, 0, result0);
    
    //printf("Before Scatter rot0: stokes x=(%e, %e, %e) y=(%e, %e, %e)", gsl_vector_get(result0_x,0), gsl_vector_get(result0_x,1), gsl_vector_get(result0_x,2), gsl_vector_get(result0_y,0), gsl_vector_get(result0_y,1), gsl_vector_get(result0_y,2));
    
    //fprintf(fPtr, "y_tilde: %e, %e, %e y_tilde_rot_result: %e, %e, %e\n", *(y_tilde+0), *(y_tilde+1), *(y_tilde+2), gsl_vector_get(y_tilde_rot_result,0), gsl_vector_get(y_tilde_rot_result,1), gsl_vector_get(y_tilde_rot_result,2));
    
    /*
     printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot0, 0,0), gsl_matrix_get(rot0, 0,1), gsl_matrix_get(rot0, 0,2));
     printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot0, 1,0), gsl_matrix_get(rot0, 1,1), gsl_matrix_get(rot0, 1,2));
     printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot0, 2,0), gsl_matrix_get(rot0, 2,1), gsl_matrix_get(rot0, 2,2));
     */
    
    //set values of ph_p_prime equal to the result and get new phi from result
    *(ph_p_prime+1)=gsl_vector_get(result0,0);
    *(ph_p_prime+2)=0;//gsl_vector_get(result,1); //just directly setting it to 0 now?
    *(ph_p_prime+3)=gsl_vector_get(result0,2);
    
    phi1=atan2(gsl_vector_get(result0,2), gsl_vector_get(result0,0));
    
    
    //printf("rotation 1: %e, %e, %e\n",  *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    //fprintf(fPtr, "Photon Phi: %e\n", phi1);
    //printf("make sure the vector view is good: %e, %e, %e,%e\n", *(ph_p_prime+0), gsl_vector_get(&ph_p.vector,0), gsl_vector_get(&ph_p.vector,1), gsl_vector_get(&ph_p.vector,2));
    
    
    //rotate around y to bring it all along x
    gsl_matrix_set(rot1, 1,1,1);
    gsl_matrix_set(rot1, 0,0,cos(-phi1));
    gsl_matrix_set(rot1, 2,2,cos(-phi1));
    gsl_matrix_set(rot1, 0,2,-sin(-phi1));
    gsl_matrix_set(rot1, 2,0,sin(-phi1));
    gsl_blas_dgemv(CblasNoTrans, 1, rot1, &ph_p.vector, 0, result1);
    
    //fprintf(fPtr, "y_tilde: %e, %e, %e y_tilde_rot vector view: %e, %e, %e\n", *(y_tilde+0), *(y_tilde+1), *(y_tilde+2), gsl_vector_get(&y_tilde_rot.vector,0), gsl_vector_get(&y_tilde_rot.vector,1), gsl_vector_get(&y_tilde_rot.vector,2));
    
    /*
     printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot1, 0,0), gsl_matrix_get(rot1, 0,1), gsl_matrix_get(rot1, 0,2));
     printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot1, 1,0), gsl_matrix_get(rot1, 1,1), gsl_matrix_get(rot1, 1,2));
     printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot1, 2,0), gsl_matrix_get(rot1, 2,1), gsl_matrix_get(rot1, 2,2));
     */
    
    //set values of ph_p_prime equal to the result and get new phi from result
    *(ph_p_prime+1)=*(ph_p_prime+0);//why setting it to the energy?
    *(ph_p_prime+2)=gsl_vector_get(result1,1);
    *(ph_p_prime+3)=0; //just directly setting it to 0 now?
    
    //printf("rotation 2: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
    
    //know that the stokes y axis is in -y_hat direction and stokes x asis is in the z_hat direction due to rotations and making inclimg photn come along x_hat direction, dont need to rotate the stokes plane/vector. this happens as the rotations occur (tested in python code)
    //double checking here
    //printf("Before Scatter: stokes x=(%e, %e, %e) y=(%e, %e, %e) ph_p=(%e, %e, %e, %e)\n", gsl_vector_get(result1_x,0), gsl_vector_get(result1_x,1), gsl_vector_get(result1_x,2), gsl_vector_get(result1_y,0), gsl_vector_get(result1_y,1), gsl_vector_get(result1_y,2), *(ph_p_prime+0), *(ph_p_prime+1), *(ph_p_prime+2), *(ph_p_prime+3));
    
    
    //determine if the scattering will occur between photon and electron
    //scattering_occured=comptonScatter(&theta, &phi, rand, fPtr); //determine the angles phi and theta for the photon to scatter into using thompson differential cross section
    scattering_occured=kleinNishinaScatter(&theta, &phi, *(ph_p_prime+0), *(s+1), *(s+2), rand, fPtr);//determine the angles phi and theta for the photon to scatter into using KN differential cross section, if the photon will end up scattering
    
    //fprintf(fPtr,"Phi: %e, Theta: %e\n", phi, theta);
    //theta=2.4475668271885342;
    //phi=4.014719957630734;
    //*(s+0)=1; //should be 1.0
    //*(s+1)=1;
    //*(s+2)=0;
    //*(s+3)=0;
    
    
    if (scattering_occured==1)
    {
        //perform scattering and compute new 4-momenta of electron and photon
        //scattered photon 4 momentum
        gsl_vector_set(result, 0, (*(ph_p_prime+0))/(1+ (( (*(ph_p_prime+0))*(1-cos(theta)) )/(M_EL*C_LIGHT )) ) ); // scattered energy of photon
        gsl_vector_set(result, 1, gsl_vector_get(result,0)*cos(theta) );
        gsl_vector_set(result, 2, gsl_vector_get(result,0)*sin(theta)*sin(phi) );//assume phi is clockwise from z to y
        gsl_vector_set(result, 3, gsl_vector_get(result,0)*sin(theta)*cos(phi) );
        //fprintf(fPtr, "New ph_p0=%e Old= %e\n", gsl_vector_get(result,0), *(ph_p_prime+0));
        //gsl_vector_fprintf(fPtr,result, "%e" );
        
        //recalc x_tilde from rotation about y by angle theta do x_tilde=y_tilde X v_ph
        //test =gsl_vector_view_array(gsl_vector_ptr(result, 1), 3);
        
        //scatt_result is a dummy, dont need to change the stokes parameters here, just need to find the axis such that y is out of the plane of k_o-k see Ito figure 12 in polarized emission from stratisfied jets
        
        //gsl_blas_ddot(&y_tilde_rot.vector, &test.vector, &dotprod_1);
        //fprintf(fPtr, "Angle between the  y_tilde_rot and the photon velocity vector is: %e\n", acos(dotprod_1/ gsl_blas_dnrm2(&test.vector))*180/M_PI);
        //gsl_vector_fprintf(fPtr,&y_tilde_rot.vector, "%e" );
        //gsl_vector_fprintf(fPtr,&x_tilde_rot.vector, "%e" );
        
        //exit(0);
        //calculate electron 4 momentum
        //prescattered photon 4 momentum
        gsl_vector_set(whole_ph_p, 0, (*(ph_p_prime+0)));
        gsl_vector_set(whole_ph_p, 1, (*(ph_p_prime+1)));
        gsl_vector_set(whole_ph_p, 2, (*(ph_p_prime+2)));
        gsl_vector_set(whole_ph_p, 3, (*(ph_p_prime+3)));
        
        gsl_vector_sub(whole_ph_p,result); //resut is saved into ph_p vector, unscattered-scattered 4 mometum of photon
        gsl_vector_add(&el_p.vector ,whole_ph_p);
        /*
         printf("After scattering:\n");
         printf("el_p: %e, %e, %e,%e\n", gsl_vector_get(&el_p.vector,0), gsl_vector_get(&el_p.vector,1), gsl_vector_get(&el_p.vector,2), gsl_vector_get(&el_p.vector,3));
         printf("ph_p: %e, %e, %e,%e\n", gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2), gsl_vector_get(result,3));
         */
        
        //rotate back to comoving frame
        *(ph_p_prime+0)=gsl_vector_get(result,0);
        *(ph_p_prime+1)=gsl_vector_get(result,1); //set values of photon prime momentum from doing the scattering to use the vector view of it in dot product
        *(ph_p_prime+2)=gsl_vector_get(result,2);
        *(ph_p_prime+3)=gsl_vector_get(result,3);
        gsl_matrix_set_all(rot1,0);
        gsl_matrix_set(rot1, 1,1,1);
        gsl_matrix_set(rot1, 0,0,cos(-phi1));
        gsl_matrix_set(rot1, 2,2,cos(-phi1));
        gsl_matrix_set(rot1, 0,2,sin(-phi1));
        gsl_matrix_set(rot1, 2,0,-sin(-phi1));
        gsl_blas_dgemv(CblasNoTrans, 1, rot1, &ph_p.vector, 0, result1);
        /*
         printf("Photon Phi: %e\n", phi1);
         printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot1, 0,0), gsl_matrix_get(rot1, 0,1), gsl_matrix_get(rot1, 0,2));
         printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot1, 1,0), gsl_matrix_get(rot1, 1,1), gsl_matrix_get(rot1, 1,2));
         printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot1, 2,0), gsl_matrix_get(rot1, 2,1), gsl_matrix_get(rot1, 2,2));
         */
        
        //set values of ph_p_prime to result1 from undoing 2nd rotation
        *(ph_p_prime+1)=gsl_vector_get(result1,0);
        *(ph_p_prime+2)=gsl_vector_get(result1,1);
        *(ph_p_prime+3)=gsl_vector_get(result1,2);
        //printf("Undo rotation 2: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
        //ignore the electron, dont care about it, undo the first rotation
        gsl_matrix_set_all(rot0,0);
        gsl_matrix_set(rot0, 2,2,1);
        gsl_matrix_set(rot0, 0,0,cos(-phi0));
        gsl_matrix_set(rot0, 1,1,cos(-phi0));
        gsl_matrix_set(rot0, 0,1,sin(-phi0));
        gsl_matrix_set(rot0, 1,0,-sin(-phi0));
        gsl_blas_dgemv(CblasNoTrans, 1, rot0, &ph_p.vector, 0, result0);
        
        
        /*
         printf("Photon Phi: %e\n", phi0);
         printf("Rotation Matrix 0: %e,%e, %e\n", gsl_matrix_get(rot0, 0,0), gsl_matrix_get(rot0, 0,1), gsl_matrix_get(rot0, 0,2));
         printf("Rotation Matrix 1: %e,%e, %e\n", gsl_matrix_get(rot0, 1,0), gsl_matrix_get(rot0, 1,1), gsl_matrix_get(rot0, 1,2));
         printf("Rotation Matrix 2: %e,%e, %e\n", gsl_matrix_get(rot0, 2,0), gsl_matrix_get(rot0, 2,1), gsl_matrix_get(rot0, 2,2));
         */
        
        //do the scattering of the stokes vector
        //rotate it by phi and then scatter it and rotate back and then renormalize it such that i=1
        //if (STOKES_SWITCH != 0)
        #if STOKES_SWITCH == ON
        {
            //orient the stokes coordinate system such that its perpendicular to the scattering plane
            findXY(gsl_vector_ptr(ph_p_orig, 1),z_axis_electron_rest_frame, x_tilde, y_tilde);
            findXY(gsl_vector_ptr(result0,0),gsl_vector_ptr(ph_p_orig, 1), x_tilde_new, y_tilde_new);
            phi=findPhi(x_tilde, y_tilde, x_tilde_new, y_tilde_new);
            mullerMatrixRotation(phi, s, fPtr);
            
            //find the theta between the incoming and scattered photons, by doing dot product and taking arccos of it
            double dot_prod_result=(gsl_vector_get(ph_p_orig,1)*gsl_vector_get(result0,0)+gsl_vector_get(ph_p_orig,2)*gsl_vector_get(result0,1)+gsl_vector_get(ph_p_orig,3)*gsl_vector_get(result0,2) )/(gsl_vector_get(ph_p_orig,0)*(*(ph_p_prime+0))) ;
            if ((dot_prod_result<-1) || (dot_prod_result>1))
            {
                //printf("The old dot poduct was %e, the new one is %e\n",dot_prod_result, round(dot_prod_result));
                dot_prod_result=round(dot_prod_result);//do this rounding so numerical error that causes value to be <-1 or >1 gets rounded and becomes a real value if its close enough to these limits
            }

            theta=acos(dot_prod_result);
            
            //do the scattering of the stokes parameters
            gsl_matrix_set(scatt, 0,0,1.0+pow(cos(theta), 2.0)+((1-cos(theta))*(gsl_vector_get(ph_p_orig,0) - gsl_vector_get(result,0))/(M_EL*C_LIGHT ) ) ); //following lundman's matrix
            gsl_matrix_set(scatt, 0,1, sin(theta)*sin(theta));
            gsl_matrix_set(scatt, 1,0, sin(theta)*sin(theta));
            gsl_matrix_set(scatt, 1,1,1.0+cos(theta)*cos(theta));
            gsl_matrix_set(scatt, 2,2, 2.0*cos(theta));
            gsl_matrix_set(scatt, 3,3, 2.0*cos(theta)+ ((cos(theta))*(1-cos(theta))*(gsl_vector_get(ph_p_orig,0) - gsl_vector_get(result,0))/(M_EL*C_LIGHT )) );
            //gsl_matrix_scale(scatt, (gsl_vector_get(result,0)/(*(ph_p_prime+0)))*((gsl_vector_get(result,0)/(*(ph_p_prime+0))))*0.5*3*THOM_X_SECT/(8*M_PI) ); //scale the matrix by 0.5*r_0^2 (\epsilon/\epsilon_0)^2 DONT NEED THIS BECAUSE WE NORMALIZE STOKES VECTOR SO THIS CANCELS ITSELF OUT
            gsl_blas_dgemv(CblasNoTrans, 1, scatt, &stokes.vector, 0, scatt_result);
            /*
             fprintf(fPtr,"before s: %e, %e, %e,%e\n", gsl_vector_get(&stokes.vector,0), gsl_vector_get(&stokes.vector,1), gsl_vector_get(&stokes.vector,2), gsl_vector_get(&stokes.vector,3));
             fprintf(fPtr,"Scatt Matrix 0: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 0,0), gsl_matrix_get(scatt, 0,1), gsl_matrix_get(scatt, 0,2), gsl_matrix_get(scatt, 0,3));
             fprintf(fPtr,"Scatt Matrix 1: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 1,0), gsl_matrix_get(scatt, 1,1), gsl_matrix_get(scatt, 1,2), gsl_matrix_get(scatt, 1,3));
             fprintf(fPtr,"Scatt Matrix 2: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 2,0), gsl_matrix_get(scatt, 2,1), gsl_matrix_get(scatt, 2,2), gsl_matrix_get(scatt, 2,3));
             fprintf(fPtr,"Scatt Matrix 3: %e,%e, %e, %e\n", gsl_matrix_get(scatt, 3,0), gsl_matrix_get(scatt, 3,1), gsl_matrix_get(scatt, 3,2), gsl_matrix_get(scatt, 3,3));
             fprintf(fPtr,"s: %e, %e, %e,%e\n", gsl_vector_get(scatt_result,0), gsl_vector_get(scatt_result,1), gsl_vector_get(scatt_result,2), gsl_vector_get(scatt_result,3));
             */
            
            
            //normalize and rotate back
            *(s+0)=gsl_vector_get(scatt_result,0)/gsl_vector_get(scatt_result,0); //should be 1.0
            *(s+1)=gsl_vector_get(scatt_result,1)/gsl_vector_get(scatt_result,0);
            *(s+2)=gsl_vector_get(scatt_result,2)/gsl_vector_get(scatt_result,0);
            *(s+3)=gsl_vector_get(scatt_result,3)/gsl_vector_get(scatt_result,0);
            //fprintf(fPtr,"s after norm: %e, %e, %e,%e\n", gsl_vector_get(&stokes.vector,0), gsl_vector_get(&stokes.vector,1), gsl_vector_get(&stokes.vector,2), gsl_vector_get(&stokes.vector,3));
            
            
            //need to find current stokes coordinate system defined in the plane of k-k_0
            findXY(gsl_vector_ptr(result0,0),gsl_vector_ptr(ph_p_orig, 1), x_tilde, y_tilde);
            
            //then find the new coordinate system between scattered photon 4 onetum and the z axis
            findXY(gsl_vector_ptr(result0,0),z_axis_electron_rest_frame, x_tilde_new, y_tilde_new);
            
            //find phi to transform between the two coodinate systems
            phi=findPhi(x_tilde, y_tilde, x_tilde_new, y_tilde_new);
            
            //do the rotation
            mullerMatrixRotation(phi, s, fPtr);
        }
        #endif
        
        //now update the array with the new scattered photon 4 monetum
        *(ph_p_prime+1)=gsl_vector_get(result0,0);
        *(ph_p_prime+2)=gsl_vector_get(result0,1);
        *(ph_p_prime+3)=gsl_vector_get(result0,2);
        
        //gsl_blas_ddot(&y_tilde_rot.vector, &ph_p.vector, &dotprod_1);
        //fprintf(fPtr, "Angle between the  y_tilde_rot and the photon velocity vector is: %e\n", acos(dotprod_1/ gsl_blas_dnrm2(&ph_p.vector))*180/M_PI);
        
        //printf("Undo rotation 1: %e, %e, %e, %e\n",  *(ph_p_prime+0), *(ph_p_prime+1),  *(ph_p_prime+2),  *(ph_p_prime+3));
        //deboost photon to lab frame
        *(negative_el_v+0)=(-1*(*(el_v+0)));
        *(negative_el_v+1)=(-1*(*(el_v+1)));
        *(negative_el_v+2)=(-1*(*(el_v+2)));
        
        lorentzBoost(negative_el_v, ph_p_prime, ph_comov, 'p', fPtr);
        //printf("Undo boost 1: %e, %e, %e, %e\n",  *(ph_comov+0), *(ph_comov+1),  *(ph_comov+2),  *(ph_comov+3));
        
        
        //dont need to find stokes vector and do previosu rotations, can just find the stokes coordinates in function because the stokes coordinate vectors rotate with the photon vector and no rotations to a new stokes coordinate system are needed
        //if (STOKES_SWITCH != 0)
        #if STOKES_SWITCH == ON
        {
            stokesRotation(negative_el_v, (ph_p_prime+1), (ph_comov+1), s, fPtr);
        }
        #endif
        
        //exit(0);
    }
    
    gsl_matrix_free(rot0); gsl_matrix_free(rot1);gsl_matrix_free(scatt);gsl_vector_free(result0);gsl_vector_free(result1);gsl_vector_free(result);
    gsl_vector_free(scatt_result);gsl_vector_free(ph_p_orig);
    gsl_vector_free(whole_ph_p);free(ph_p_prime);free(el_p_prime);free(el_v); free(negative_el_v); free(z_axis_electron_rest_frame);
    
    return scattering_occured;
}

int comptonScatter(double *theta, double *phi, gsl_rng * rand, FILE *fPtr)
{
    
        double y_dum, f_x_dum, x_dum;
        
        //generate random theta and phi angles for scattering
        *phi=gsl_rng_uniform(rand)*2*M_PI;
        //printf("Phi: %e\n", phi);
    
        y_dum=1; //initalize loop to get a random theta
        f_x_dum=0;
        while (y_dum>f_x_dum)
        {
            y_dum=gsl_rng_uniform(rand)*1.09;
            x_dum=gsl_rng_uniform(rand)*M_PI;
            f_x_dum=sin(x_dum)*(1+cos(x_dum)*cos(x_dum));
        }
        *theta=x_dum;
        
        return 1;
}

int kleinNishinaScatter(double *theta, double *phi, double p0, double q, double u, gsl_rng * rand, FILE *fPtr)
{
    //sample theta using:  https://doi.org/10.13182/NSE11-57
    double phi_dum=0, cos_theta_dum=0, f_phi_dum=0, f_cos_theta_dum=0, f_theta_dum=0, phi_y_dum=0, cos_theta_y_dum=0, KN_x_section_over_thomson_x_section=0, rand_num=0;
    double mu=0, phi_max=0, norm=0;
    int will_scatter=0;
    double energy_ratio=  p0/(M_EL*C_LIGHT ); //h*nu / mc^2 , units of p0 is erg/c
    
    //determine the KN cross section over the thomson cross section
    KN_x_section_over_thomson_x_section= kleinNishinaCrossSection(energy_ratio);
    rand_num=gsl_rng_uniform(rand);
        
    if (rand_num<= KN_x_section_over_thomson_x_section)
    {
        //fprintf(fPtr,"In If!\n");
        //fflush(fPtr);
    
        //sample a theta and phi from the differential cross sections
        phi_y_dum=1; //initalize loop to get a random phi and theta
        cos_theta_y_dum=1;
        f_cos_theta_dum=0;
        f_phi_dum=0;
        
        while ((cos_theta_y_dum>f_cos_theta_dum))
        {
            //do phi and theta seperately, sample theta using:  https://doi.org/10.13182/NSE11-57
            cos_theta_y_dum=gsl_rng_uniform(rand)*2;
            cos_theta_dum=gsl_rng_uniform(rand)*2-1;
            f_cos_theta_dum=pow((1+energy_ratio*(1-cos_theta_dum)),-2)*(energy_ratio*(1-cos_theta_dum)+(1/(1+energy_ratio*(1-cos_theta_dum))) + cos_theta_dum*cos_theta_dum);
            
        }
        *theta=acos(cos_theta_dum);
        mu=1+energy_ratio*(1-cos(*theta));
        f_theta_dum=(pow(mu, -1.0) + pow(mu, -3.0) - pow(mu, -2.0)*sin(*theta)*sin(*theta))*sin(*theta);
        
        while ((phi_y_dum>f_phi_dum) )
        {
            
            #if STOKES_SWITCH == OFF
            {
                //not considering polarization therefore can jjst sample between 0 and 2*pi evenly
                phi_dum=gsl_rng_uniform(rand)*2*M_PI;
                phi_y_dum=-1; // this is to exit the while statement
                
                //fprintf(fPtr," phi_dum: %e\n", phi_dum);
                //fflush(fPtr);

            }
            #else
            {
                if (u==0 && q==0)
                {
                    phi_dum=gsl_rng_uniform(rand)*2*M_PI;
                    phi_y_dum=-1; // this is to exit the while statement

                }
                else
                {
                    //if we are considering polarization calulate the norm for the distributiion to be between 1 and 0
                    phi_max=fabs(atan2(-u,q))/2.0;
                    norm=(f_theta_dum + pow(mu, -2.0)*sin(*theta)*sin(*theta)*sin(*theta) * (q*cos(2*phi_max)-u*sin(2*phi_max)));
                    //fprintf(fPtr,"norm: %e\n", norm);
                    //fflush(fPtr);
                    
                    phi_y_dum=gsl_rng_uniform(rand);
                    phi_dum=gsl_rng_uniform(rand)*2*M_PI;
                    f_phi_dum=(f_theta_dum + pow(mu, -2.0)*sin(*theta)*sin(*theta)*sin(*theta) * (q*cos(2*phi_dum)-u*sin(2*phi_dum)))/norm; //signs on q and u based on Lundman/ McMaster
                    
                    //fprintf(fPtr,"phi_y_dum: %e, theta_dum: %e, mu: %e, f_theta_dum: %e, phi_dum: %e, f_phi_dum: %e, u: %e, q: %e\n", phi_y_dum, theta_dum, mu, f_theta_dum, phi_dum, f_phi_dum, u, q);
                    //fflush(fPtr);

                }
            }
            #endif
            
        }
        *phi=phi_dum;
        
        will_scatter=1;
    }
    else
    {
        will_scatter=0;
    }
    
    return will_scatter;
}

double kleinNishinaCrossSection(double energy_ratio)
{
    /*
        Calculate the total cross section normalized by the thompson cross section
        given the photon energy (in the electron rest frame) normalized by the electron rest mass energy.
        This is given in the grmonty code and the low energy limit is more continuous and gives very small error
        with respect to the full cross section value (needed for integration when calculating modified cross section).
        Verified that this form of the KN cross section reproduces Rybicki and lightman KN cross section formula
    */
    double result=0;

    if (energy_ratio >= 1e-3)
    {
        result=(3. / 4.) * (2. / (energy_ratio * energy_ratio) +
                 (1. / (2. * energy_ratio) -
                  (1. + energy_ratio) / (energy_ratio * energy_ratio * energy_ratio)) * log(1. +
                                2. * energy_ratio) +
                 (1. + energy_ratio) / ((1. + 2. * energy_ratio) * (1. + 2. * energy_ratio))
        );
    }
    else
    {
        result=(1. - 2. * energy_ratio);
    }

    return result;
}

double scatteredPhotonWeight(double weight, double bias, double optical_depth)
{
    #if NONTHERMAL_E_DIST != OFF
        if (bias != 1)
        {
            return weight*(-gsl_expm1(-optical_depth))/(-gsl_expm1(-bias*optical_depth));
            //return weight*(1-exp(-optical_depth))/(1-exp(-bias*optical_depth));
        }
        else
        {
            return weight;
        }
    #else
        return weight;
    #endif
}
