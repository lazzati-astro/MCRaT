//
// Created by Tyler Parsotan on 11/8/25.
//

#include "mcrat.h"

double calculateOpticalDepth(struct photon *ph, struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    int ph_block_index;
    double tau=0;
    double ph_phi=0;
    double fl_v_x=0, fl_v_y=0, fl_v_z=0; //to hold the fluid velocity in MCRaT coordinates
    double ph_v_norm=0, fl_v_norm=0;
    double n_cosangle=0, n_dens_lab_tmp=0;
    double beta=0, fluid_beta[3];


    ph_block_index=ph->nearest_block_index;

    //save values
    n_dens_lab_tmp = (hydro_data->dens_lab)[ph_block_index];

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

    tau = M_P/(n_dens_lab_tmp))/THOM_X_SECT/(1.0-beta*n_cosangle);

    return tau;
}