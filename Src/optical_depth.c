//
// Created by Tyler Parsotan on 11/8/25.
//

#include "mcrat.h"

double calculateOpticalDepth(struct photon *ph, struct hydro_dataframe *hydro_data, gsl_rng *rand, FILE *fPtr)
{
    int ph_block_index, i;
    double tau=0;
    double ph_phi=0;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates
    double ph_v_norm=0, fl_v_norm=0;
    double n_cosangle=0, thermal_n_dens_lab=0;
    double beta=0, fluid_beta[3], fluid_factor=0;
    #if NONTHERMAL_E_DIST == OFF
        double norm_cross_section=0;
    #else
        double norm_cross_section[1+N_GAMMA];
    #endif


    ph_block_index=ph->nearest_block_index;

    //save values
    thermal_n_dens_lab = (hydro_data->dens_lab)[ph_block_index]/M_P;

    #if DIMENSIONS == THREE
        hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[ph_block_index], (hydro_data->v1)[ph_block_index], (hydro_data->v2)[ph_block_index], (hydro_data->r0)[ph_block_index], (hydro_data->r1)[ph_block_index], (hydro_data->r2)[ph_block_index]);
    #elif DIMENSIONS == TWO_POINT_FIVE
        ph_phi=atan2(ph->r1, ph->r0);
        hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[ph_block_index], (hydro_data->v1)[ph_block_index], (hydro_data->v2)[ph_block_index], (hydro_data->r0)[ph_block_index], (hydro_data->r1)[ph_block_index], ph_phi);
    #else
        ph_phi=atan2(ph->r1, ph->r0);
        //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
        hydroVectorToCartesian(&fluid_beta, (hydro_data->v0)[ph_block_index], (hydro_data->v1)[ph_block_index], 0, (hydro_data->r0)[ph_block_index], (hydro_data->r1)[ph_block_index], ph_phi);
    #endif

    fl_v_x = fluid_beta[0];
    fl_v_y = fluid_beta[1];
    fl_v_z = fluid_beta[2];

    fl_v_norm = sqrt(fl_v_x*fl_v_x+fl_v_y*fl_v_y+fl_v_z*fl_v_z);
    ph_v_norm = sqrt((ph->p1)*(ph->p1)+(ph->p2)*(ph->p2)+(ph->p3)*(ph->p3));

    //(*(n_cosangle+i))=((fl_v_x* (ph->p1))+(fl_v_y* (ph->p2))+(fl_v_z* (ph->p3)))/(fl_v_norm*ph_v_norm ); //find cosine of the angle between the photon and the fluid velocities via a dot product
    n_cosangle = ((fl_v_x* (ph->p1))+(fl_v_y* (ph->p2))+(fl_v_z* (ph->p3)))/(fl_v_norm*ph_v_norm ); //make 1 for cylindrical otherwise its undefined

    beta = sqrt(1.0-1.0/((hydro_data->gamma)[ph_block_index]*(hydro_data->gamma)[ph_block_index]));

    fluid_factor=(1.0-beta*n_cosangle);

    //TODO: extend this to the non-thermal electron dist
    #if NONTHERMAL_E_DIST == OFF
        getCrossSection( ph->comv_p0,  (hydro_data->temp)[ph_block_index], &norm_cross_section,  rand, fPtr);
        tau = 1/(thermal_n_dens_lab)/(THOM_X_SECT*norm_cross_section)/fluid_factor;
    #else
        getCrossSection( ph->comv_p0,  (hydro_data->temp)[ph_block_index], norm_cross_section,  rand, fPtr);

        //calculate the thermal tau
        tau = 1/(thermal_n_dens_lab)/(THOM_X_SECT*(*(norm_cross_section+0))/fluid_factor;

        //calculate the nonthermal tau
    #endif

    return tau;
}

void getCrossSection(double photon_comv_e, double fluid_temp, double *cross_section, gsl_rng *rand, FILE *fPtr)
{
    //this returns the cross section normalized by the thompson cross section
    #if TAU_CALCULATION == TABLE
        *(cross_section+0)=getThermalCrossSection(photon_comv_e, fluid_temp, rand, fPtr);
        #if NONTHERMAL_E_DIST != OFF
            getNonThermalCrossSection( photon_comv_e, (cross_section+1), rand, fPtr);
        #endif
    #else
        //if we are directly calcualting the optical depth, just use the thompson cross section
        *(cross_section+0)=1;
    #endif

}

double getThermalCrossSection(double photon_comv_e, double fluid_temp, gsl_rng *rand, FILE *fPtr)
{
    //this returns the thermal cross section normalized by the thompson cross section

    double result=0;

    #if TAU_CALCULATION == TABLE
        double normalized_photon_comv_e=photon_comv_e/(M_EL*C_LIGHT ); //h*nu / mc^2 , units of p0 is erg/c
        double theta=fluid_temp*(K_B/(M_EL*C_LIGHT*C_LIGHT ));
        //fprintf(fPtr, "normalized_photon_comv_e: %e, theta: %e\n", normalized_photon_comv_e, theta);
        //fflush(fPtr);
        result = pow(10.0, interpolateThermalHotCrossSection(log10(normalized_photon_comv_e), log10(theta), rand, fPtr));
    #else
        result = 1;
    #endif

    return result;
}

#if NONTHERMAL_E_DIST != OFF
    double getNonThermalCrossSection(double photon_comv_e, double *subgroup_interpolated_results, gsl_rng *rand, FILE *fPtr)
    {
        double normalized_photon_comv_e=photon_comv_e/(M_EL*C_LIGHT ); //h*nu / mc^2 , units of p0 is erg/c

        interpolateSubgroupNonThermalHotCrossSection(log10(normalized_photon_comv_e), subgroup_interpolated_results, rand, fPtr);
        fprintf(fPtr, "NonThermal test: %g %g %g %g\n", log10(normalized_photon_comv_e), *(subgroup_interpolated_results+0), *(subgroup_interpolated_results+1), *(subgroup_interpolated_results+2));

    }
#endif

