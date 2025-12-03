void initalizePhotonList(struct photonList *photon_list);

void freePhotonList(struct photonList *photon_list);

void allocatePhotonListMemory(struct photonList *photon_list, int n_photons);

void reallocatePhotonListMemory(struct photonList *photon_list, int new_capacity);

void appendToPhotonList(struct photonList *photon_list, struct photon *ph);
