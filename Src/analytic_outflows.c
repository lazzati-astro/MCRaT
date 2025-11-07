//
// Created by Tyler Parsotan on 11/7/25.
//

#include "mcrat.h"

void cylindricalPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
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

void sphericalPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    double  gamma_infinity=100, lumi=1e54, r00=1e8; //shopuld be 10^57
    //double  gamma_infinity=5, lumi=1e52, r00=1e8; //shopuld be 10^57
    double vel=0, r=0;
    int i=0;

    fprintf(fPtr, "The Spherical Outflow values are: Gamma_infinity=%e, Luminosity=%e erg/s, r_0=%e cm \n", gamma_infinity, lumi, r00);
    fflush(fPtr);

    for (i=0; i<hydro_data->num_elements; i++)
    {
        if (((hydro_data->r))[i] >= (r00*gamma_infinity))
        {
            ((hydro_data->gamma))[i]=gamma_infinity;
            ((hydro_data->pres))[i]=(lumi*pow(r00, 2.0/3.0)*pow(((hydro_data->r))[i], -8.0/3.0) )/(12.0*M_PI*C_LIGHT*pow(gamma_infinity, 4.0/3.0));
        }
        else
        {
            ((hydro_data->gamma))[i]=((hydro_data->r))[i]/r00;
            ((hydro_data->pres))[i]=(lumi*pow(r00, 2.0))/(12.0*M_PI*C_LIGHT*pow(((hydro_data->r))[i], 4.0) );
        }

        ((hydro_data->dens))[i]=lumi/(4*M_PI*pow(((hydro_data->r))[i], 2.0)*pow(C_LIGHT, 3.0)*gamma_infinity*(((hydro_data->gamma))[i]));
        ((hydro_data->dens_lab))[i]=(((hydro_data->dens))[i])*(((hydro_data->gamma))[i]);
        ((hydro_data->temp))[i]=pow(3*(((hydro_data->pres))[i])/(A_RAD) ,1.0/4.0);

        vel=sqrt(1-(pow(((hydro_data->gamma))[i], -2.0)));

        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE

            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2));
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
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2)+pow(((hydro_data->r2))[i], 2));
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
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r2))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r; //need to figure this out
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif

        #endif

        //fprintf(fPtr,"Gamma: %lf\nR: %lf\nPres: %e\nvel %lf\nX: %lf\nY %lf\nVx: %lf\nVy: %lf\nDens: %e\nLab_Dens: %e\nTemp: %lf\n", *(gamma+i), *(r+i), *(pres+i), vel, *(x+i), *(y+i), *(vx+i), *(vy+i), *(dens+i), *(dens_lab+i), *(temp+i));
    }

}

void structuredFireballPrep(struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    //This model is provided by Lundman, Peer, Ryde 2014, use this to compare our MCRaT polarization to their polarizations
    double  gamma_0=100, lumi=1e52, r00=1e8, theta_j=1e-2, p=4; //theta_j in paper is 1e-2, 3e-2, 1e-1 and p is 1,2,4
    double T_0=pow(lumi/(4*M_PI*r00*r00*A_RAD*C_LIGHT), 1.0/4.0);
    double eta=0, r_sat=0, r;
    double vel=0, theta_ratio=0;
    int i=0;

    fprintf(fPtr, "The Structured Spherical Outflow values are: Gamma_0=%e, Luminosity=%e erg/s, r_0=%e cm, theta_j=%e rad, p=%e \n", gamma_0, lumi, r00, theta_j, p);
    fflush(fPtr);

    for (i=0; i<hydro_data->num_elements; i++)
    {


        theta_ratio=((hydro_data->theta)[i])/theta_j;
        eta=gamma_0/sqrt(1+pow(theta_ratio, 2*p));

        if ((hydro_data->theta)[i] >= theta_j*pow(gamma_0/2, 1.0/p))
        {
            //*(gamma+i)=2; //outside with of shear layer have gamma be 2 like in paper
            eta=2.0;
        }

        r_sat=eta*r00;

        if (((hydro_data->r)[i]) >= r_sat)
        {
            (hydro_data->gamma)[i]=eta;
            (hydro_data->temp)[i]=T_0*pow(r_sat/((hydro_data->r)[i]), 2.0/3.0)/eta;
        }
        else
        {
            (hydro_data->gamma)[i]=((hydro_data->r)[i])/r_sat; //not sure if this is right but it shouldn't matter since we're injecting our photons far from r00
            (hydro_data->temp)[i]=T_0;
        }

        vel=sqrt(1-(pow((hydro_data->gamma)[i], -2.0)));
        (hydro_data->dens)[i] = M_P*lumi/(4*M_PI*M_P*C_LIGHT*C_LIGHT*C_LIGHT*eta*vel*((hydro_data->gamma)[i])*((hydro_data->r)[i])*((hydro_data->r)[i])); //equation paper has extra c, but then units dont work out
        (hydro_data->dens_lab)[i]=((hydro_data->dens)[i])*((hydro_data->gamma)[i]);
        (hydro_data->pres)[i]=(A_RAD*pow((hydro_data->temp)[i], 4.0))/(3);

        #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE

            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2));
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
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r1))[i], 2)+pow(((hydro_data->r2))[i], 2));
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
                r=sqrt(pow(((hydro_data->r0))[i], 2)+ pow(((hydro_data->r2))[i], 2));
                ((hydro_data->v0))[i]=(vel*(((hydro_data->r0))[i]))/r;
                ((hydro_data->v1))[i]=0;
                ((hydro_data->v2))[i]=(vel*(((hydro_data->r2))[i]))/r;
            #endif

        #endif

        //fprintf(fPtr,"eta: %lf\nr_sat: %lf\nGamma: %lf\nR: %lf\nTheta: %lf\nPres: %e\nvel %lf\nX: %lf\nY %lf\nVx: %lf\nVy: %lf\nDens: %e\nLab_Dens: %e\nTemp: %lf\n\n", eta, r_sat, *(gamma+i), *(r+i), (*(theta+i)), *(pres+i), vel, *(x+i), *(y+i), *(vx+i), *(vy+i), *(dens+i), *(dens_lab+i), *(temp+i));

    }

}
