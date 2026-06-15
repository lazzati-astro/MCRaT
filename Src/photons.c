#include "mcrat.h"

void initalizePhotonList(struct photonList *photon_list)
{
    //initialize pointers in photon_list to NULL for debugging
    photon_list->photons=NULL;
    photon_list->sorted_indexes=NULL;
    
    //initalize the number of photons, number of null photons, and the list capacity to 0
    photon_list->num_photons=0;
    photon_list->num_null_photons=0;
    photon_list->list_capacity=0;

}

void freePhotonList(struct photonList *photon_list)
{
    free(photon_list->photons);
    free(photon_list->sorted_indexes);
    photon_list->photons=NULL;
    photon_list->sorted_indexes=NULL;
    photon_list->num_photons=0;
    photon_list->list_capacity=0;
    photon_list->num_null_photons=0;

}

void allocatePhotonListMemory(struct photonList *photon_list, int n_photons)
{
    photon_list->photons = malloc (n_photons * sizeof (struct photon ));
    photon_list->sorted_indexes = malloc(n_photons*sizeof(int));
    
    photon_list->list_capacity=n_photons;
    
}

void reallocatePhotonListMemory(struct photonList *photon_list, int new_capacity)
{
    //extend the photon list to be new_capacity elements long
    int i=0, old_list_capacity=photon_list->list_capacity;
    struct photon *new_photons = realloc(photon_list->photons, new_capacity * sizeof(struct photon));
    int *new_sorted_indexes = realloc(photon_list->sorted_indexes, new_capacity * sizeof(int));
    
    //make sure that the realloc worked
    if (new_photons != NULL)
    {
        /* everything ok */
        photon_list->photons = new_photons;
    }
    else
    {
        /* problems!!!! */
        printf("Error with reserving space to hold new photons\n");
        exit(1);
    }

    if (new_sorted_indexes != NULL)
    {
        /* everything ok */
        photon_list->sorted_indexes = new_sorted_indexes;
    }
    else
    {
        /* problems!!!! */
        printf("Error with reserving space to hold new sorted index array\n");
        exit(1);
    }
    photon_list->list_capacity=new_capacity;
    
    /*
     we only do this when there are no more null photons to replace with real ones and we need to grow the list. so we can artificially say that all the newly added photons are real and in the next loop we make them all null photons explicity and keep track of them appropriately
     */
    photon_list->num_photons+=(new_capacity-old_list_capacity);
    
    //assign all new elements to be null photons, setNullPhoton decreases num_photons by 1 each time it is called
    for (i=old_list_capacity;i<new_capacity;i++)
    {
        setNullPhoton(photon_list, i);
    }
}

void setPhotonList(struct photonList *photon_list, struct photon *ph_array, int num_photons)
{
    int i=0, null_photon_count=0;
    //this copies an array of photons into the Photon list struct. This overwrites any prior list that was saved in the struct.
    if (photon_list->photons != NULL)
    {
        freePhotonList(photon_list);
    }
    allocatePhotonListMemory(photon_list, num_photons);
    
    memcpy(photon_list->photons, ph_array, num_photons*sizeof(struct photon));
    photon_list->list_capacity=num_photons;
    photon_list->num_photons=num_photons;
    photon_list->num_null_photons=0;
    
    //get the actual number of null photons in the array we dont want to assume that there are 0.
    for (i = 0; i < num_photons; i++)
    {
        if (photon_list->photons[i].type == NULL_PHOTON)
        {
            incrementNullPhotonNum(photon_list);
        }
    }

}

void addToPhotonList(struct photonList *photon_list, struct photon *ph, size_t num_photons)
{
    int idx=0, i=0, j=0, new_capacity=0;
    int *null_ph_indexes=NULL;
    
    //add a photon to the photonList photons array
    // If list is full, and we have no null photons to fill in then double capacity
    if ((photon_list->num_photons >= photon_list->list_capacity) || (photon_list->num_null_photons <= num_photons))
    {
        if (photon_list->list_capacity * 2 > photon_list->list_capacity + num_photons)
        {
            new_capacity = photon_list->list_capacity * 2;
        }
        else
        {
            //take care of the case where doubling the memory is not enough
            new_capacity = photon_list->list_capacity * (double) ceil((double) (photon_list->list_capacity + num_photons)/ (double) photon_list->list_capacity) ;
        }
        reallocatePhotonListMemory(photon_list, new_capacity);
    }
    
    //if we have just 1 photons to append we can just assign it to end of list or first location of NULL  photon
    if (num_photons == 1)
    {
        //if we have no null photons, just append the photon to the list
        if (photon_list->num_null_photons == 0)
        {
            idx=photon_list->num_photons;
            
        }
        else
        {
            //we need to find the null photon index and overwrite the photon there
            for (i = 0; i < photon_list->list_capacity; i++)
            {
                if (photon_list->photons[i].type == NULL_PHOTON)
                {
                    idx=i;
                    i=photon_list->list_capacity;
                }
            }
        }
        
        // Copy photon into list
        memcpy(&photon_list->photons[idx], ph, sizeof(struct photon));
        incrementPhotonNum(photon_list);

        //photon_list->num_photons++;
        //photon_list->num_null_photons--;

    }
    else
    {
        //if we have many phtons we are appending, we want to find any/all null photons and assign all the photons to these spots. Note when we reallocate memory above, we set all new photons to be null
        null_ph_indexes=malloc((photon_list->num_null_photons)*sizeof(int));
        
        //verify that num of null phtons is greater than number of photons to append
        if  (num_photons>(photon_list->num_null_photons))
        {
            printf("Adding to the photon list has failed. the number of null photons in the list is less than the number of photons to add to the list. %d vs %d", (photon_list->num_null_photons), num_photons );
            exit(1);
        }
        //TODO: this finds all the null photons even if we only need some number that is less than photon_list->num_null_photons, this can be optimized.
        for (i=0;i<photon_list->list_capacity; i++)
        {
            if (photon_list->photons[i].type == NULL_PHOTON)
            {
                // if the weight is 0, this is a photons that has been absorbed and is now null
                *(null_ph_indexes+j)=i;
                j++;
                //fprintf(fPtr, "NULL PHOTON INDEX %d\n", i);
                //fflush(fPtr);
                
            }
        }
        
        //now we have the null photon indexes, we can assign the photons to these indexes
        for (i=0;i<num_photons;i++)
        {
            //only do the assigning for non-null photons
            if ((ph+i)->type != NULL_PHOTON)
            {
                idx=(*(null_ph_indexes+i));
                // Copy photon into list
                memcpy(&photon_list->photons[idx], (ph+i), sizeof(struct photon));
                incrementPhotonNum(photon_list);
                //photon_list->num_photons++;
                //photon_list->num_null_photons--;
            }
        }
        
        free(null_ph_indexes);

        
    }
    
}

void setNullPhoton(struct photonList *photon_list, int index)
{
    int i=0;
    
    //we set a photon in the photon list at index to be NULL
    photon_list->photons[index].type = NULL_PHOTON;
    photon_list->photons[index].weight = 0;
    photon_list->photons[index].nearest_block_index = -1;
    photon_list->photons[index].recalc_properties = 0;
    
    //valgrind wants us to initalize the other members of the struct
    photon_list->photons[index].p0=0;
    photon_list->photons[index].p1=0;
    photon_list->photons[index].p2=0;
    photon_list->photons[index].p3=0;
    photon_list->photons[index].comv_p0=0;
    photon_list->photons[index].comv_p1=0;
    photon_list->photons[index].comv_p2=0;
    photon_list->photons[index].comv_p3=0;
    photon_list->photons[index].r0=0;
    photon_list->photons[index].r1=0;
    photon_list->photons[index].r2=0;
    photon_list->photons[index].s0=0;
    photon_list->photons[index].s1=0;
    photon_list->photons[index].s2=0;
    photon_list->photons[index].s3=0;
    photon_list->photons[index].num_scatt=0;
    photon_list->photons[index].total_scattering_opacity=0;
    photon_list->photons[index].time_to_scatter = FLT_MAX / C_LIGHT;
    
    #if NONTHERMAL_E_DIST != OFF
        for (i=0;i<(1+N_GAMMA);i++)
        {
            photon_list->photons[index].scattering_opacity[i]=0; //the optical depths that are calculated for thermal + non-thermal electrons with the nonthermal electron subgroups
            photon_list->photons[index].scattering_bias[i]=0;
        }
    #endif
    
    
    #if SYNCHROTRON_SWITCH == ON
        photon_list->photons[index].absorption_opacity=0;
    #endif
    
    //photon_list->num_null_photons++;
    //photon_list->num_photons--;
    incrementNullPhotonNum(photon_list);

}


struct photon* getPhoton(struct photonList *photon_list, int index)
{
    return &photon_list->photons[index];
}

void incrementPhotonNum(struct photonList *photon_list)
{
    (photon_list->num_photons)+=1;
    (photon_list->num_null_photons)-=1;

    verifyPhotonNum(photon_list);

}

void incrementNullPhotonNum(struct photonList *photon_list)
{
    (photon_list->num_photons)-=1;
    (photon_list->num_null_photons)+=1;
    
    verifyPhotonNum(photon_list);
    
}

void verifyPhotonNum(struct photonList *photon_list)
{
    if (photon_list->num_photons+photon_list->num_null_photons != photon_list->list_capacity)
    {
        printf("Error with incrementing real (%d) or null (%d) photons in the photonList of capacity (%d). Conservation of photon error\n", photon_list->num_photons, photon_list->num_null_photons, photon_list->list_capacity);
        exit(1);

    }
}

struct photon createPhoton()
{
    //this function is something that the user can call in creating their custom photon injection algorithm. The photons are explicitly not set to be a null-type photon (see setNullPhoton above) since that is actually used in the code algorithm, instead we set all the fields to be NAN as is applicable so we can easily test the fields that the user has over written and load then into the photon that will be saved into the photonList struct that gets used throughout the code.
    struct photon ph;
    initializePhoton(&ph);
    return ph;
}


void initializePhoton(struct photon *ph)
{
    if (ph == NULL) return;

    ph->type = '\0';

    ph->p0 = NAN;
    ph->p1 = NAN;
    ph->p2 = NAN;
    ph->p3 = NAN;

    ph->comv_p0 = NAN;
    ph->comv_p1 = NAN;
    ph->comv_p2 = NAN;
    ph->comv_p3 = NAN;

    ph->r0 = NAN;
    ph->r1 = NAN;
    ph->r2 = NAN;

    ph->s0 = NAN;
    ph->s1 = NAN;
    ph->s2 = NAN;
    ph->s3 = NAN;

    ph->num_scatt = NAN;
    ph->recalc_properties = -1;
    ph->weight = NAN;
    ph->nearest_block_index = -1;
    ph->time_to_scatter = NAN;

    #if SCATTERING_BIAS_SWITCH != OFF
        for (int i = 0; i < 1 + N_GAMMA; ++i)
        {
            ph->scattering_opacity[i] = NAN;
            ph->scattering_bias[i] = NAN;
        }
    #endif

    ph->total_scattering_opacity = NAN;
    
    #if SYNCHROTRON_SWITCH == ON
        ph->absorption_opacity=NAN; //this holds the opacity for synchrotron absorption (multiply by path length to get optical depth when modifying the weight correspondent with synchrotron self absorption
    #endif

}

double samplePhotonPhi(gsl_rng * rand, FILE *fPtr)
{
    return gsl_rng_uniform(rand)*2*M_PI;
}

double samplePhotonTheta(double *velocity, gsl_rng * rand, FILE *fPtr)
{
    //trying to overwrite com_v_theta based on sampling of lab anisotropic angle distribution of photons
    //see eg Section 3.2.1 @ doi.org/10.1088/0004-637X/807/1/31 & Section 3.5 @ doi.org/10.3847/1538-4357/ac75cb
    // and section 6.2 here: Nordin Nobuoka, J. 2025, SPIRO: a code that couples Monte Carlo photons to relativistic hydrodynamics - Applications to hot astrophysical plasmas, https://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-368279
    double com_v_theta;
    gsl_vector_view b=gsl_vector_view_array(velocity, 3);
    double beta=gsl_blas_dnrm2(&b.vector);
    double y_dum=1; //initalize loop
    double yfr_dum=0;
    while (y_dum>yfr_dum)
    {
        com_v_theta=2*gsl_rng_uniform_pos(rand)-1; //cos(angle) is from -1 to 1
        //printf("%lf, %lf ",gsl_rng_uniform_pos(rand), (*(temps+i)));
        y_dum=gsl_rng_uniform_pos(rand);
         
        yfr_dum=0.5*(1+beta*com_v_theta); //propability density of angle of photon with respect to fluid motion (doppler boosting factor)
    }
    com_v_theta=acos(com_v_theta);
    
    return com_v_theta;

}

void saveUserDefinePhoton(struct photon *ph_orig, struct photon *ph_user, struct hydro_dataframe *hydro_data, int hydro_index, gsl_rng * rand, FILE *fPtr)
{
    //here we go through each member of the photon struct and see if the user has filled it in with non-NAN value. If so, assign it to the original photon which will be saved to the photonList that is used in MCRaT
    //curently we dont let the photon weight, the position, num_scatt, nearest_block_index, recalc_properties, photon_type, and optical depth related stuff be changed by the user
    
    double orig_energy = 0, position_phi=NAN;
    double *p = NULL;
    double *boost = NULL;
    double *l_boost = NULL;
    bool recalc_lab_momentum = false, recalc_fluid_momentum = false;
    
    //right now we dont allow the user to redefine the number of scatterings, the weights, the position, the nearest_block_index, the recalc_properties, or any optical depth properties. If they do try to do that print out an error and exit
    if (ph_user->type != '\0')
    {
        fprintf(fPtr, "The user cannot redefine the injected photon type.\n");
        fflush(fPtr);
        exit(1);
    }
    
    if (!isnan(ph_user->weight))
    {
        fprintf(fPtr, "The user cannot redefine the injected photon weight.\n");
        fflush(fPtr);
        exit(1);
    }
    
    if (!isnan(ph_user->num_scatt))
    {
        fprintf(fPtr, "The user cannot redefine the injected photon number of scatterings.\n");
        fflush(fPtr);
        exit(1);
    }

    if (ph_user->nearest_block_index != -1)
    {
        fprintf(fPtr, "The user cannot redefine the injected photon's injected block index field.\n");
        fflush(fPtr);
        exit(1);

    }


    if (ph_user->recalc_properties != -1)
    {
        fprintf(fPtr, "The user cannot redefine the injected photon's recalc properties switch.\n");
        fflush(fPtr);
        exit(1);

    }
    
    if (!isnan(ph_user->total_scattering_opacity))
    {
        fprintf(fPtr, "The user cannot redefine the injected photon's total optical depth field.\n");
        fflush(fPtr);
        exit(1);

    }

    if (!isnan(ph_user->time_to_scatter))
    {
        fprintf(fPtr, "The user cannot redefine the injected photon's time to scatter.\n");
        fflush(fPtr);
        exit(1);

    }
    
    if (!isnan(ph_user->r0) || !isnan(ph_user->r1) || !isnan(ph_user->r2))
    {
        fprintf(fPtr, "The user cannot redefine the injected photon's position.\n");
        fflush(fPtr);
        exit(1);

    }
        
    #if SCATTERING_BIAS_SWITCH != OFF
        for (int i = 0; i < 1 + N_GAMMA; ++i)
        {
            if (!isnan(ph_user->scattering_opacity[i]))
            {
                fprintf(fPtr, "The user cannot redefine the injected photon's optical depth array.\n");
                fflush(fPtr);
                exit(1);

            }
            if (!isnan(ph_user->scattering_bias[i]))
            {
                fprintf(fPtr, "The user cannot redefine the injected photon's scattering bias switch.\n");
                fflush(fPtr);
                exit(1);

            }

        }
    #endif


    p = malloc(4*sizeof(double));
    boost = malloc(3*sizeof(double));
    l_boost = malloc(4*sizeof(double));

    
    //here we allow the user to just set an energy in the comoving frame and we overwrite the orig 4 momentum with that energy if the other 3 elements of the ph_user comoving 4 momentum are nans
    if (!isnan(ph_user->comv_p0) && isnan(ph_user->comv_p1) && isnan(ph_user->comv_p2) && isnan(ph_user->comv_p3))
    {
        orig_energy = ph_orig->comv_p0;
        ph_orig->comv_p0 = ph_user->comv_p0;
        ph_orig->comv_p1 *= (ph_user->comv_p0/orig_energy);
        ph_orig->comv_p2 *= (ph_user->comv_p0/orig_energy);
        ph_orig->comv_p3 *= (ph_user->comv_p0/orig_energy);
        recalc_lab_momentum=true;
    }
        
    //see if the user completely overwrites the comoving 4 momentum
    if (!isnan(ph_user->comv_p0) && !isnan(ph_user->comv_p1) && !isnan(ph_user->comv_p2) && !isnan(ph_user->comv_p3))
    {
        ph_orig->comv_p0 = ph_user->comv_p0;
        ph_orig->comv_p1 = ph_user->comv_p1;
        ph_orig->comv_p2 = ph_user->comv_p2;
        ph_orig->comv_p3 = ph_user->comv_p3;
        recalc_lab_momentum=true;

    }

        
    //here we allow the user to just set an energy in the lab frame and we overwrite the orig 4 momentum with that energy if the other 3 elements of the ph_user lab 4 momentum are nans
    if (!isnan(ph_user->p0) && isnan(ph_user->p1) && isnan(ph_user->p2) && isnan(ph_user->p3))
    {
        orig_energy = ph_orig->p0;
        ph_orig->p0 = ph_user->p0;
        ph_orig->p1 *= (ph_user->p0/orig_energy);
        ph_orig->p2 *= (ph_user->p0/orig_energy);
        ph_orig->p3 *= (ph_user->p0/orig_energy);
        recalc_fluid_momentum=true;

    }
        
    //see if the user completely overwrites the lab 4 momentum
    if (!isnan(ph_user->p0) && !isnan(ph_user->p1) && !isnan(ph_user->p2) && !isnan(ph_user->p3))
    {
        ph_orig->p0 = ph_user->p0;
        ph_orig->p1 = ph_user->p1;
        ph_orig->p2 = ph_user->p2;
        ph_orig->p3 = ph_user->p3;
        recalc_fluid_momentum=true;

    }
        
    //see if the user wants to overwrite the stokes, remember this is stokes in the lab frame
    if (!isnan(ph_user->s0) && !isnan(ph_user->s1) && !isnan(ph_user->s2) && !isnan(ph_user->s3))
    {
        //check if the s0 ==1 and that the other 3 added in quadrature up to one
        if ((gsl_fcmp(ph_user->s0, 1.0, GSL_DBL_EPSILON) == 0) && (ph_user->s0*ph_user->s0 >= ph_user->s1*ph_user->s1 + ph_user->s2*ph_user->s2 + ph_user->s3*ph_user->s3))
        {
            ph_orig->s0 = ph_user->s0;
            ph_orig->s1 = ph_user->s1;
            ph_orig->s2 = ph_user->s2;
            ph_orig->s3 = ph_user->s3;
        }
        else
        {
            fprintf(fPtr, "The user defined stokes parameters are not the expected input for MCRaT. S0 should be 1 and S1^2+S2^2+S3^2>=S1^2. user input: S0: %lf, S1: %lf, S2: %lf, S3: %lf\n", ph_user->s0, ph_user->s1, ph_user->s2, ph_user->s3);
            fflush(fPtr);
        }

    }
    
    //since we arent modifying the position of the photons based on the user generated photon, we still use the original photon position to get it's phi
    position_phi=fmod(atan2(ph_orig->r1, ph_orig->r0) * 180.0/M_PI + 360.0, 360.0) *M_PI/180;
    
    if (recalc_lab_momentum)
    {
        //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
        #if DIMENSIONS == THREE
            hydroVectorToCartesian(boost, (hydro_data->v0)[hydro_index], (hydro_data->v1)[hydro_index], (hydro_data->v2)[hydro_index], (hydro_data->r0)[hydro_index], (hydro_data->r1)[hydro_index], (hydro_data->r2)[hydro_index]);
        #elif DIMENSIONS == TWO_POINT_FIVE
            hydroVectorToCartesian(boost, (hydro_data->v0)[hydro_index], (hydro_data->v1)[hydro_index], (hydro_data->v2)[hydro_index], (hydro_data->r0)[hydro_index], (hydro_data->r1)[hydro_index], position_phi);
        #else
            //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
            hydroVectorToCartesian(boost, (hydro_data->v0)[hydro_index], (hydro_data->v1)[hydro_index], 0, (hydro_data->r0)[hydro_index], (hydro_data->r1)[hydro_index], position_phi);
        #endif
        (*(boost+0))*=-1;
        (*(boost+1))*=-1;
        (*(boost+2))*=-1;
        
        *(p+0)=ph_orig->comv_p0;
        *(p+1)=ph_orig->comv_p1;
        *(p+2)=ph_orig->comv_p2;
        *(p+3)=ph_orig->comv_p3;
            
        //boost to lab frame
        lorentzBoost(boost, p, l_boost, 'p', fPtr);

        ph_orig->p0=(*(l_boost+0));
        ph_orig->p1=(*(l_boost+1));
        ph_orig->p2=(*(l_boost+2));
        ph_orig->p3=(*(l_boost+3));

    }

    if (recalc_fluid_momentum)
    {
        //populate boost matrix, not sure why multiplying by -1, seems to give correct answer in old python code...
        #if DIMENSIONS == THREE
            hydroVectorToCartesian(boost, (hydro_data->v0)[hydro_index], (hydro_data->v1)[hydro_index], (hydro_data->v2)[hydro_index], (hydro_data->r0)[hydro_index], (hydro_data->r1)[hydro_index], (hydro_data->r2)[hydro_index]);
        #elif DIMENSIONS == TWO_POINT_FIVE
            hydroVectorToCartesian(boost, (hydro_data->v0)[hydro_index], (hydro_data->v1)[hydro_index], (hydro_data->v2)[hydro_index], (hydro_data->r0)[hydro_index], (hydro_data->r1)[hydro_index], position_phi);
        #else
            //this may have to change if PLUTO can save vectors in 3D when conidering 2D sim
            hydroVectorToCartesian(boost, (hydro_data->v0)[hydro_index], (hydro_data->v1)[hydro_index], 0, (hydro_data->r0)[hydro_index], (hydro_data->r1)[hydro_index], position_phi);
        #endif
        
        *(p+0)=ph_orig->comv_p0;
        *(p+1)=ph_orig->comv_p1;
        *(p+2)=ph_orig->comv_p2;
        *(p+3)=ph_orig->comv_p3;

        
        //boost to fluid frame
        lorentzBoost(boost, p, l_boost, 'p', fPtr);

        ph_orig->comv_p0=(*(l_boost+0));
        ph_orig->comv_p1=(*(l_boost+1));
        ph_orig->comv_p2=(*(l_boost+2));
        ph_orig->comv_p3=(*(l_boost+3));

        
    }
        
    free(p);free(boost);free(l_boost);

    
    
}
