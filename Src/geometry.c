//
//  geometry.c
//  
//  functions that convert between different coordinates, calculates different volume elements, etc. Anything related to the geometry/coordinates of the hydro simulation vs. MCRaT cartesian coordinates
//  Created by Tyler Parsotan on 7/23/21.
//

#include "mcrat.h"

void mcratCoordinateToHydroCoordinate(double *ph_hydro_coord, double mcrat_r0, double mcrat_r1, double mcrat_r2)
{
    //function to convert MCRaT cartesian coordinate in 3D to the proper hydro coordinates
    double r0=-1, r1=-1, r2=-1;
    
    #if DIMENSIONS == 2
        //can be cartesian (x,z), polar (r,phi), spherical (r,theta), cylindrical (r,z),
        //for cartesian have to reduce MCRaT x,y to hydro x coord and leave z
        //(for polar have to reduce MCRaT x,y,z to radius and calculate phi from projected x axis), NOT SUPPORTING POLAR 2D
        //for spherical do the same as for polar
        //for cylindrical do the same as cartesian
        //ASSUME AXISYMMETRIC JET AXIS ALONG SECOND COORDINATE
    
        #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
            r0=sqrt(mcrat_r0*mcrat_r0+mcrat_r1*mcrat_r1); //x coordinate or r coordinate
            r1=mcrat_r2; //coordinate along jet axis
        #endif
    
        #if GEOMETRY == SPHERICAL
            r0=sqrt(mcrat_r0*mcrat_r0+mcrat_r1*mcrat_r1+mcrat_r2*mcrat_r2); //r coordinate
            r1= acos(mcrat_r2/r0); //cooridinate along jet axis
        #endif
        
    #else
        //can be cartesian (x,y,z;same as MCRaT), spherical (r, theta, phi), Polar (r, phi, z)
        #if GEOMETRY == CARTESIAN
            r0=mcrat_r0; //x coordinate
            r1=mcrat_r1; // y coordinate
            r2=mcrat_r2; // z coordinate
        #endif
    
        #if GEOMETRY == SPHERICAL
            r0=sqrt(mcrat_r0*mcrat_r0+mcrat_r1*mcrat_r1+mcrat_r2*mcrat_r2); // r coordinate
            r1= acos(mcrat_r2/r0); //theta cooridinate along jet axis
            r2=atan2(mcrat_r1, mcrat_r0); // phi coordinate
        #endif

        #if GEOMETRY == POLAR
            r0=sqrt(mcrat_r0*mcrat_r0+mcrat_r1*mcrat_r1); // r coordinate
            r1=atan2(mcrat_r1, mcrat_r0); //phi coordinate along jet axis
            r2=mcrat_r2; // z coordinate
        #endif

    #endif
    
    *(ph_hydro_coord+0)=r0;
    *(ph_hydro_coord+1)=r1;
    *(ph_hydro_coord+2)=r2;
    
}

void hydroCoordinateToSpherical(double *r, double *theta, double r0, double r1, double r2)
{
    //this function converts hydro coordinates to spherical r and theta coordinates
    int i=0;
    double sph_r=0, sph_theta=0;//sph_theta is measured from the assumed jet axis
    
    #if DIMENSIONS == 2
        
        #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
            sph_r=pow( r0*r0+r1*r1, 0.5);
            sph_theta=atan2( r0 , r1 );

        #endif

        #if GEOMETRY == SPHERICAL
            sph_r=r0;
            sph_theta=r1;
        #endif
        
    #else

        #if GEOMETRY == CARTESIAN
            sph_r=pow(  r0 * r0 + r1 * r1 + r2 * r2 , 0.5);
            sph_theta=acos(  r2 /sph_r );
        #endif

        #if GEOMETRY == SPHERICAL
            sph_r=r0;
            sph_theta=r1;
        #endif

        #if GEOMETRY == POLAR
            sph_r=pow(  r0 * r0 + r2 * r2 , 0.5);//sqrt(r^2+z^2)
            sph_theta=acos( r2/sph_r );
        #endif

    #endif
    
    *r=sph_r;
    *theta=sph_theta;
    
}

void hydroCoordinateToMcratCoordinate(double *hydro_mcrat_coord, double hydro_r0, double hydro_r1, double hydro_r2)
{
    //converts hydro coordinate to MCRaT coordinate system (3D cartesian)
    // in 2D for cylindrical, cartesian, and spherical coordinates can just pass in the photon phi for the value of hydro_r2 since we are assuming axisymmetry here anyways
    
    double x=0, y=0, z=0;
    
    #if DIMENSIONS == 2
        
        #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
            x=hydro_r0*cos(hydro_r2);
            y=hydro_r0*sin(hydro_r2);
            z=hydro_r1;
        #endif

        #if GEOMETRY == SPHERICAL
            x=hydro_r0*sin(hydro_r1)*cos(hydro_r2);
            y=hydro_r0*sin(hydro_r1)*sin(hydro_r2);
            z=hydro_r0*cos(hydro_r1);
        #endif
        
    #else

        #if GEOMETRY == CARTESIAN
            x=hydro_r0;
            y=hydro_r1;
            z=hydro_r2;
        #endif

        #if GEOMETRY == SPHERICAL
            x=hydro_r0*sin(hydro_r1)*cos(hydro_r2);
            y=hydro_r0*sin(hydro_r1)*sin(hydro_r2);
            z=hydro_r0*cos(hydro_r1);
        #endif

        #if GEOMETRY == POLAR
            x=hydro_r0*cos(hydro_r1);
            y=hydro_r0*sin(hydro_r1);
            z=hydro_r2;
        #endif

    #endif

    *(hydro_mcrat_coord+0)=x;
    *(hydro_mcrat_coord+1)=y;
    *(hydro_mcrat_coord+2)=z;
}

void fillHydroCoordinateToSpherical(struct hydro_dataframe *hydro_data)
{
    //this function fills in the r and theta values in the hydro_data struct, which is used for photon injection, overwriting values, etc
    int i=0;
    double sph_r=0, sph_theta=0;//sph_theta is measured from the assumed jet axis
    
    for (i=0;i<hydro_data->num_elements;i++)
    {
        #if DIMENSIONS == 3
            hydroCoordinateToSpherical(&sph_r, &sph_theta, ((hydro_data->r0))[i], ((hydro_data->r1))[i], ((hydro_data->r2))[i]);
        #else
            hydroCoordinateToSpherical(&sph_r, &sph_theta, ((hydro_data->r0))[i], ((hydro_data->r1))[i], 0);
        #endif
        ((hydro_data->r))[i]=sph_r;
        ((hydro_data->theta))[i]=sph_theta;

    }
    
}

void hydroVectorToCartesian(double *cartesian_vector_3d, double v0, double v1, double v2, double x0, double x1, double x2)
{
    //takes the vector <v0, v1, v2> at position (x0, x1, x2) in the hydro coordinate system and converts it to a 3D vector
    // in 2D for cylindrical, cartesian, and spherical coordinates can just pass in the photon phi for the value of x2 since we are assuming axisymmetry here anyways
    //may need to modify if PLUTO allows for saving full 3D vectors even when it only considers 2D sims
    double transformed_vector0=0, transformed_vector1=0, transformed_vector2=0;
    
    #if DIMENSIONS == 2
        
        #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
            transformed_vector0=v0*cos(x2); //x coordinate
            transformed_vector1=v0*sin(x2) ; //y
            transformed_vector2=v1; //z
        #endif

        #if GEOMETRY == SPHERICAL
            v2=0;//no phi hat component of the vector in 2D
            transformed_vector0=v0*sin(x1)*cos(x2)+v1*cos(x1)*cos(x2)-v2*sin(x2); //x coordinate
            transformed_vector1=v0*sin(x1)*sin(x2)+v1*cos(x1)*sin(x2)+v2*cos(x2); //y
            transformed_vector2=v0*cos(x1)-v1*sin(x1); //z
        #endif
        
    #else
    
        #if GEOMETRY == CARTESIAN
            transformed_vector0=v0; //x coordinate
            transformed_vector1=v1 ; //y
            transformed_vector2=v2; //z
        #endif

        #if GEOMETRY == SPHERICAL
            transformed_vector0=v0*sin(x1)*cos(x2)+v1*cos(x1)*cos(x2)-v2*sin(x2); //x coordinate
            transformed_vector1=v0*sin(x1)*sin(x2)+v1*cos(x1)*sin(x2)+v2*cos(x2); //y
            transformed_vector2=v0*cos(x1)-v1*sin(x1); //z
        #endif
    
        #if GEOMETRY == POLAR
            transformed_vector0=v0*cos(x1)-v1*sin(x1); //x coordinate
            transformed_vector1=v0*sin(x1)+v1*cos(x1); //y
            transformed_vector2=v2; //z
        #endif

    #endif
    
    *(cartesian_vector_3d+0)=transformed_vector0;
    *(cartesian_vector_3d+1)=transformed_vector1;
    *(cartesian_vector_3d+2)=transformed_vector2;
    
}

double hydroElementVolume(struct hydro_dataframe *hydro_data, int index)
{
    //calculate the volume of a hydro element assuming axissymmetry in 2D case
    double V=0, r0_min=0, r0_max=0, r1_min=0, r1_max=0, r2_min=0, r2_max=0;
    
    r0_max=(hydro_data->r0)[index]+0.5*(hydro_data->r0_size)[index];
    r0_min=(hydro_data->r0)[index]-0.5*(hydro_data->r0_size)[index];
    r1_max=(hydro_data->r1)[index]+0.5*(hydro_data->r1_size)[index];
    r1_min=(hydro_data->r1)[index]-0.5*(hydro_data->r1_size)[index];

    #if DIMENSIONS == 2
        
        #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
            V=M_PI*(r0_max*r0_max-r0_min*r0_min)*(hydro_data->r1_size)[index];
        #endif

        #if GEOMETRY == SPHERICAL
            V=(2.0*M_PI/3.0)*(r0_max*r0_max*r0_max-r0_min*r0_min*r0_min)*(cos(r1_min)-cos(r1_max));//dV=r^2sin(theta)dr dtheta dphi
        #endif
        
    #else
    
        r2_max=(hydro_data->r2)[index]+0.5*(hydro_data->r2_size)[index];
        r2_min=(hydro_data->r2)[index]-0.5*(hydro_data->r2_size)[index];

        #if GEOMETRY == CARTESIAN
            V=(hydro_data->r0_size)[index]*(hydro_data->r1_size)[index]*(hydro_data->r2_size)[index];
        #endif

        #if GEOMETRY == SPHERICAL
            V=(1.0/3.0)*(r0_max*r0_max*r0_max-r0_min*r0_min*r0_min)*(cos(r1_min)-cos(r1_max))*(r2_max-r2_min);//dV=r^2sin(theta)dr dtheta dphi
        #endif

        #if GEOMETRY == POLAR
            V=0.5*(r0_max*r0_max-r0_min*r0_min)*(hydro_data->r1_size)[index]*(hydro_data->r2_size)[index];//dV=rdr dphi dz
        #endif

    #endif

    
    return V;
}

int findNearestBlock(int array_num, double ph_x, double ph_y, double ph_z, double *x, double  *y, double *z)
{
    double dist=0, dist_min=1e15, block_dist=0;
    int min_index=0, j=0;
    
    dist_min=1e15;//set dist to impossible value to make sure at least first distance calulated is saved
    block_dist=3e9;
    while (dist_min==1e15) //if this is true, then the algorithm hasnt found blocks within the acceptable range given by block_dist
    {
        
        for(j=0;j<array_num;j++)
        {
            //if the distance between them is within 3e9, to restrict number of possible calculations,  calulate the total distance between the box and photon
            #if DIMENSIONS == 2
                if ((fabs(ph_x- (*(x+j)))<block_dist) && (fabs(ph_y- (*(y+j)))<block_dist))
                {
                    
                    dist= pow(pow(ph_x- (*(x+j)), 2.0) + pow(ph_y- (*(y+j)) , 2.0),0.5);
                    //fprintf(fPtr,"Dist calculated as: %e, index: %d\n", dist, j);
                    //printf("In outer if statement, OLD: %e, %d\n", dist_min, min_index);
                    
                    if((dist<dist_min))
                    {
                        //fprintf(fPtr,"In innermost if statement, OLD: %e, %d\n", dist_min, min_index);
                        dist_min=dist; //save new minimum distance
                        min_index=j; //save index
                        //printf("New Min dist: %e, New min Index: %d, Array_Num: %d\n", dist_min, min_index, array_num);
                        
                    }
                    
                }
            #elif DIMENSIONS == 3
                if ((fabs(ph_x- (*(x+j)))<block_dist) && (fabs(ph_y- (*(y+j)))<block_dist) && (fabs(ph_z- (*(z+j)))<block_dist))
                {
                    dist= pow(pow(ph_x- (*(x+j)), 2.0) + pow(ph_y- (*(y+j)),2.0 ) + pow(ph_z- (*(z+j)) , 2.0),0.5);
                    if((dist<dist_min))
                    {
                        //printf("In innermost if statement, OLD: %e, %d\n", dist_min, min_index);
                        dist_min=dist; //save new minimum distance
                        min_index=j; //save index
                        //fprintf(fPtr,"New Min dist: %e, New min Index: %d, Array_Num: %e\n", dist_min, min_index, array_num);
                    }
                }
            #endif
        }
        block_dist*=10; //increase size of accepted distances for gris points, if dist_min==1e12 then the next time the acceptance range wil be larger
        
    }
    
    return min_index;
}

int findContainingBlock(double ph_hydro_r0, double ph_hydro_r1, double ph_hydro_r2, struct hydro_dataframe *hydro_data, FILE *fPtr)
{
    int i=0, within_block_index=0;
    bool is_in_block=0; //boolean to determine if the photon is outside of a grid
    
    //can parallelize here to save time?
    for (i=0;i<hydro_data->num_elements;i++)
    {
        
        is_in_block=checkInBlock(ph_hydro_r0, ph_hydro_r1, ph_hydro_r2, hydro_data, i); //(i,  ph_x,  ph_y,  ph_z,  x,   y, z,  szx,  szy);
        
            if (is_in_block)
            {
                within_block_index=i;
                //change for loop index once the block is found so the code doesnt search the rest of the grids to see if the photon is within those grids
                i=hydro_data->num_elements;
            }
        
    }
    //printf("Within Block Index:  %d\n",within_block_index);
    #if SIM_SWITCH == RIKEN || DIMENSIONS == 3
    {
        fprintf(fPtr, "3D switch is: %d and SIM switch is: %d\n", DIMENSIONS, SIM_SWITCH);
    }
    #endif
    
    if (is_in_block==0)
    {
        #if DIMENSIONS == 2
            fprintf(fPtr, "MCRaT Couldn't find a block for the photon located at r0=%e r1=%e\n", ph_hydro_r0, ph_hydro_r1);
        #else
            fprintf(fPtr, "MCRaT Couldn't find a block for the photon located at r0=%e r1=%e r2=%e in the hydro simulation coordinate system.\n", ph_hydro_r0, ph_hydro_r1, ph_hydro_r2);
        #endif
        fflush(fPtr);
        within_block_index=-1;
    }
    
    return within_block_index;
}


int checkInBlock(double ph_hydro_r0, double ph_hydro_r1, double ph_hydro_r2, struct hydro_dataframe *hydro_data, int block_index)
{
    bool is_in_block=0; //boolean to determine if the photon is outside of its previously noted block
    double x0=0, x1=0, x2=0, sz_x0=0, sz_x1=0, sz_x2=0; //coordinate and sizes of grid block, in cartesian its x,y,z in spherical its r,theta,phi
    int return_val=0;

    
    #if DIMENSIONS == 2
        is_in_block= (2*fabs( ph_hydro_r0 - (hydro_data->r0)[block_index]) - (hydro_data->r0_size)[block_index] <= 0) && (2*fabs(ph_hydro_r1 - (hydro_data->r1)[block_index] ) - (hydro_data->r1_size)[block_index]  <= 0);
    #else
        is_in_block= (2*fabs( ph_hydro_r0 - (hydro_data->r0)[block_index]) - (hydro_data->r0_size)[block_index] <= 0) && (2*fabs(ph_hydro_r1 - (hydro_data->r1)[block_index] ) - (hydro_data->r1_size)[block_index]  <= 0) && (2*fabs(ph_hydro_r2 - (hydro_data->r2)[block_index] ) - (hydro_data->r2_size)[block_index]  <= 0);
    #endif
        /*
        else
        {
            if (riken_switch==1)
            {
                
                x0=pow(pow((*(x+block_index)), 2.0) + pow((*(y+block_index)),2.0 ) + pow((*(z+block_index)) , 2.0),0.5);
                x1=acos((*(z+block_index))/pow(pow((*(x+block_index)), 2.0) + pow((*(y+block_index)),2.0 ) + pow((*(z+block_index)) , 2.0),0.5));
                x2=atan2((*(y+block_index)), (*(x+block_index)));
                
                sz_x0=(*(szy+block_index));
                sz_x1=(*(szx+block_index));
                sz_x2=(*(szx+block_index));
                
                is_in_block= (fabs(pow(pow( ph_x, 2.0) + pow(ph_y, 2.0)+pow(ph_z, 2.0),0.5) - x0) <= sz_x0/2.0) &&  (fabs(acos(ph_z/pow(pow(ph_x, 2.0) + pow(ph_y,2.0 ) + pow(ph_z , 2.0),0.5)) - x1 ) <= sz_x1/2.0)  && (fabs(atan2(ph_y, ph_x) - x2 ) <= sz_x2/2.0);
                //not sure why the code was going to this line above here for spherical test
                 
             }
        }
        */
        
    if (is_in_block)
    {
        return_val=1;
    }
    else
    {
        return_val=0;
    }
    
    return return_val;
}
