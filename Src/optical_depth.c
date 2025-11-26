//
// Created by Tyler Parsotan on 11/8/25.
//

#include "mcrat.h"

void calculateOpticalDepth(struct photon *ph, struct hydro_dataframe *hydro_data, gsl_rng *rand, FILE *fPtr)
{
    int ph_block_index, i;
    double tau=0; //this returns the total optical depth (includes thermal and nonthermal photons)
    double ph_phi=0;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates
    double ph_v_norm=0, fl_v_norm=0;
    double n_cosangle=0, thermal_n_dens_lab=0;
    double beta=0, fluid_beta[3], fluid_factor=0;
    #if NONTHERMAL_E_DIST == OFF
        double norm_cross_section=0;
    #else
        double norm_cross_section[1+N_GAMMA], nonthermal_n_dens_lab_i, nonthermal_n_dens_lab;
    #endif


    ph_block_index=ph->nearest_block_index;


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

    //save values
    thermal_n_dens_lab = (hydro_data->dens_lab)[ph_block_index]/M_P;

    //TODO: extend this to the non-thermal electron dist
    #if NONTHERMAL_E_DIST == OFF
        getCrossSection( ph->comv_p0,  (hydro_data->temp)[ph_block_index], &norm_cross_section,  rand, fPtr);
        (ph->total_optical_depth) = (thermal_n_dens_lab)*(THOM_X_SECT*norm_cross_section)*fluid_factor;
        tau = (ph->total_optical_depth);
    #else
        getCrossSection( ph->comv_p0,  (hydro_data->temp)[ph_block_index], norm_cross_section,  rand, fPtr);

        //calculate the thermal tau
        (ph->optical_depths)[0] = (thermal_n_dens_lab)*(THOM_X_SECT*(*(norm_cross_section+0)))*fluid_factor;
        fprintf(fPtr, "thermal tau: %e\n", (ph->optical_depths)[0] );

        //get the nonthermal electron density based on magnetic energy density and electron distribution
        //then multiply by gamma to get the nonthermal electron density in lab frame
        nonthermal_n_dens_lab=(hydro_data->nonthermal_dens)[ph_block_index]*(hydro_data->gamma)[ph_block_index]; // calculateNonthermalElectronDens(hydro_data, ph_block_index)*(hydro_data->gamma)[ph_block_index] ;

        //calculate the nonthermal tau
        for (i=0;i<N_GAMMA;i++)
        {
            nonthermal_n_dens_lab_i=nonthermal_n_dens_lab*(hydro_data->electron_dens_subgroup)[i];
            (ph->optical_depths)[i+1] = 1/(nonthermal_n_dens_lab_i)/(THOM_X_SECT*(*(norm_cross_section+(i+1))))/fluid_factor;
            fprintf(fPtr, "nonthermal_n_dens_lab_i: %e, subgroup_dens: %e, norm_cross_section: %e ith tau: %e\n", nonthermal_n_dens_lab_i, (hydro_data->electron_dens_subgroup)[i], *(norm_cross_section+(i+1)), (ph->optical_depths)[i+1] );
            fflush(fPtr);

        }

        tau=0;
        for (i=0;i<N_GAMMA+1;i++)
        {
            tau += (ph->optical_depths)[i];
        }
        (ph->total_optical_depth) = tau;

        fprintf(fPtr, "total tau: %e\n", (ph->total_optical_depth) );
        fflush(fPtr);

    #endif

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
        double theta=calcDimlessTheta(fluid_temp);
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
        int i=0;
        double normalized_photon_comv_e=photon_comv_e/(M_EL*C_LIGHT ); //h*nu / mc^2 , units of p0 is erg/c

        interpolateSubgroupNonThermalHotCrossSection(log10(normalized_photon_comv_e), subgroup_interpolated_results, rand, fPtr);
        //fprintf(fPtr, "NonThermal test: %g %g %g %g\n", log10(normalized_photon_comv_e), *(subgroup_interpolated_results+0), *(subgroup_interpolated_results+1), *(subgroup_interpolated_results+2));

        //the interpolated results are returned as log(cross section/thompson cross section) so we want to raise to
        //power 10 to get back to the ratio  cross section/thompson cross section that we expect from this function
        for (i=0;i<N_GAMMA+1;i++)
        {
            *(subgroup_interpolated_results+i) = pow(10, *(subgroup_interpolated_results+i));
        }

    }
#endif

double calculateThermalScatteringBias(double alpha_parameter, double average_dimless_theta, double cell_dimless_theta, double tau)
{
    return fmax(1.0, alpha_parameter*cell_dimless_theta/(average_dimless_theta*tau));
}

double calculateNonthermalScatteringBias(double thermal_scatt_bias, double thermal_tau, double nonthermal_tau)
{
    return thermal_scatt_bias*thermal_tau/nonthermal_tau;
}
