#include "mcrat.h"

void initalizePhotonList(struct photonList *photon_list)
{
    //initialize pointers in photon_list to NULL for debugging
    photon_list->photons=NULL;
    
    //initalize the number of photons, number of null photons, and the list capacity to 0
    photon_list->num_photons=0;
    photon_list->num_null_photons=0;
    photon_list->list_capacity=0;

}

void freePhotonList(struct photonList *photon_list)
{
    free(photon_list->photons);
    photon_list->photons=NULL;
}

void allocatePhotonListMemory(struct photonList *photon_list, int n_photons)
{
    photon_list->photons=malloc (n_photons * sizeof (struct photon ));
    
    photon_list->list_capacity=n_photons;
    
}

void reallocatePhotonListMemory(struct photonList *photon_list, int new_capacity)
{
    //extend the photon list to be new_capacity elements long
    struct photon *new_photons = realloc(photon_list->photons, new_capacity * sizeof(struct photon));
    
    //make sure that the realloc worked
    if (new_photons != NULL)
    {
        /* everything ok */
        photon_list->photons = new_photons;
        photon_list->list_capacity=new_capacity;
    }
    else
    {
        /* problems!!!! */
        printf("Error with reserving space to hold new photons\n");
        exit(1);
    }
    
}

void addToPhotonList(struct photonList *photon_list, struct photon *ph)
{
    int idx=0, i=0;
    
    //add a photon to the photonList photons array
    // If list is full, and we have no null photons to fill in then double capacity
    if ((photon_list->num_photons >= photon_list->list_capacity) && (photon_list->num_null_photons == 0))
    {
        int new_capacity = photon_list->list_capacity * 2;
        reallocatePhotonListMemory(photon_list, new_capacity);
    }
    
    //if we have no null photons, just append the photon to the list
    if (photon_list->num_null_photons == 0)
    {
        idx=photon_list->num_photons;
        
    }
    else
    {
        //we need to find the null photon index and overwrite the photon there
        for (i = 0; i < photon_list->num_photons; i++)
        {
            if (photon_list->photons[i].type == NULL_PHOTON)
            {
                idx=i;
                i=photon_list->num_photons;
            }
        }
        photon_list->num_null_photons--;
    }
    
    // Copy photon into list
    memcpy(&(photon_list->photons)[idx], ph, sizeof(struct photon));
    photon_list->num_photons++;
    
}




