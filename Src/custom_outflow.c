//
//  custom_outflow.c
//  
//
//  Created by Tyler Parsotan on 7/14/22.
//

#include "mcrat.h"

void customOutflowPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    double  gamma_infinity=100, t_comov=1e5, ddensity=3e-7;// the comoving temperature in Kelvin, and the comoving density in g/cm^2
    int i=0;
    double vel=sqrt(1-pow(gamma_infinity, -2.0)), lab_dens=gamma_infinity*ddensity;
    
    fprintf(fPtr, "The Cylindrical Outflow values are: Gamma_infinity=%e, T_comv=%e K, comv dens=%e g/cm^3 \n", gamma_infinity, t_comov, ddensity);
    fflush(fPtr);
    
    for (i=0; i<hydro_data->num_elements; i++)
    {
        ((hydro_data->gamma))[i]=gamma_infinity;
        ((hydro_data->dens))[i]=ddensity;
        ((hydro_data->dens_lab))[i]=lab_dens;
        ((hydro_data->pres))[i]=(A_RAD*pow(t_comov, 4.0))/(3);
        ((hydro_data->temp))[i]=t_comov; //just assign t_comov
        
        
        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
            
            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                ((hydro_data->v0))[i]=0;
                ((hydro_data->v1))[i]=vel; //geometry dependent want this to be parallel to jet axis
            #endif

            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel*cos(((hydro_data->r1))[i]);//rhat
                ((hydro_data->v1))[i]=-vel*sin(((hydro_data->r1))[i]);//theta hat direction
            #endif
        
            #if DIMENSIONS == TWO_POINT_FIVE
                //have to make sure that the 3rd vctro direction is set to 0 in 2.5D case
                ((hydro_data->v2))[i]=0;
            #endif
            
        #else

            #if GEOMETRY == CARTESIAN
                ((hydro_data->v0))[i]=0;
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=vel;
            #endif


            #if GEOMETRY == SPHERICAL
                ((hydro_data->v0))[i]=vel*cos(((hydro_data->r1))[i]);//rhat
                ((hydro_data->v1))[i]=-vel*sin(((hydro_data->r1))[i]);//theta hat direction
                ((hydro_data->v2))[i]=0;
            #endif

            #if GEOMETRY == POLAR
                ((hydro_data->v0))[i]=0;
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=vel;
            #endif

        #endif
        

    }
    
}

