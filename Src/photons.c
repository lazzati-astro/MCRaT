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
    
    //get the actual number of null photons in the array we dont want to assume that there are 0.
    for (i = 0; i < num_photons; i++)
    {
        if (photon_list->photons[i].type == NULL_PHOTON)
        {
            null_photon_count++;
        }
    }

    
    photon_list->num_null_photons=null_photon_count;
    
}

void addToPhotonList(struct photonList *photon_list, struct photon *ph, size_t num_photons)
{
    int idx=0, i=0, j=0, new_capacity=0;
    int *null_ph_indexes=NULL;
    
    //add a photon to the photonList photons array
    // If list is full, and we have no null photons to fill in then double capacity
    if ((photon_list->num_photons >= photon_list->list_capacity) && (photon_list->num_null_photons <= num_photons))
    {
        if (photon_list->list_capacity * 2 > photon_list->list_capacity + num_photons)
        {
            new_capacity = photon_list->list_capacity * 2;
        }
        else
        {
            //take care of the case where doubling the memory is not enough
            new_capacity = photon_list->list_capacity * (num_photons/photon_list->list_capacity) ;
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
    photon_list->photons[index].total_optical_depth=0;
    
    #if NONTHERMAL_E_DIST != OFF
    for (i=0;i<(1+N_GAMMA);i++)
    {
        photon_list->photons[index].optical_depths[i]=0; //the optical depths that are calculated for thermal + non-thermal electrons with the nonthermal electron subgroups
        photon_list->photons[index].scattering_bias[i]=0;
    }
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
        printf("Error with incremenitng real or null photon in the photonList. conservation of photon error\n");
        exit(1);

    }
}
