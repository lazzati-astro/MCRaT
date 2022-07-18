//
//  custom_outflow.c
//  
//
//  Created by Tyler Parsotan on 7/14/22.
//

#include "mcrat.h"

void customOutflowPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    //for now the velocity is in units of c, but will need to change this everywhere once I am done testing
    double  beta=2.0, t_comov=3*1.60218e-9/K_B, m_dot_out=0.8, tau_0=2;// the comoving temperature of 3keV in Kelvin, and the initial unitless m_dot, initial optical depth
    double r_0=1e9; //choose some radius at base of outflow in cm, also is inner limit of simulations
    double N_0=tau_0*(beta-1)/(THOM_X_SECT*r_0); //initial number density of particles
    double vel=0, rho=0, r=0, gamma;
    int i=0;
    
    fprintf(fPtr, "The Custom Outflow values are: beta=%e, T_comv=%e K, m_dot=%e, tau_0=%e, r_0=%e \n", beta, t_comov, m_dot_out, tau_0, r_0);
    fflush(fPtr);
    
    for (i=0; i<hydro_data->num_elements; i++)
    {
        r=((hydro_data->r))[i];
        vel=m_dot_out*pow(r_0/r, 2-beta); //C_LIGHT*, see first comment about velocity in units of c
        
        //should this be the lab frame density or the fluid frame density?
        //this may not matter since its a factor of ~2
        rho=N_0*M_P*pow(r_0/r, beta); 
        gamma=1/sqrt(1-pow(vel/C_LIGHT, 2.0));
        
        
        ((hydro_data->gamma))[i]=gamma;
        ((hydro_data->dens))[i]=rho;
        ((hydro_data->dens_lab))[i]=rho*gamma;
        ((hydro_data->pres))[i]=(A_RAD*pow(t_comov, 4.0))/(3);
        ((hydro_data->temp))[i]=t_comov; //just assign t_comov
        
        
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
            
            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=(vel*(((hydro_data->r1))[i]))/r; //geometry dependent want this to be radial
            #endif

            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel;//rhat
                ((hydro_data->v1))[i]=0;//theta hat direction
            #endif

            #if DIMENSIONS == TWO_POINT_FIVE
                //have to make sure that the 3rd vctro direction is set to 0 in 2.5D case
                ((hydro_data->v2))[i]=0;
            #endif
            
        #else

            #if GEOMETRY == CARTESIAN
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=(vel*(((hydro_data->r1))[i]))/r; //geometry dependent want this to be radial
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif


            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel;//rhat
                ((hydro_data->v1))[i]=0;//theta hat direction
                ((hydro_data->v2))[i]=0;
            #endif

            #if GEOMETRY == POLAR
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r; //need to figure this out
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif

        #endif
        

    }
    
}

