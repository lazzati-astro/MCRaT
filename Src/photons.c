#include "mcrat.h"

void initalizePhotonList(struct photonList *photon_list)
{
    //initialize pointers in photon_list to NULL for debugging
    photon_list->photons=NULL;
    
    //initalize the number of photons and the lsit capacity to 0
    photon_list->num_photons=0;
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

void appendToPhotonList(struct photonList *photon_list, struct photon *ph)
{
    //add a photon to the photonList photons array
    // If list is full, double capacity
    if (photon_list->num_photons >= photon_list->list_capacity)
    {
        int new_capacity = photon_list->list_capacity * 2;
        photon_list_resize(photon_list, new_capacity);
    }
    
    // Copy photon into list
    memcpy(&photon_list->photons[plist->num_photons], ph, sizeof(struct photon));
    photon_list->num_photons++;
    
}


