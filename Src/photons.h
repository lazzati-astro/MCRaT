void initalizePhotonList(struct photonList *photon_list);

void freePhotonList(struct photonList *photon_list);

void allocatePhotonListMemory(struct photonList *photon_list, int n_photons);

void reallocatePhotonListMemory(struct photonList *photon_list, int new_capacity);

void setPhotonList(struct photonList *photon_list, struct photon *ph_array, int num_photons);

void addToPhotonList(struct photonList *photon_list, struct photon *ph, size_t num_photons);

void setNullPhoton(struct photonList *photon_list, int index);

struct photon* getPhoton(struct photonList *photon_list, int index);

void incrementPhotonNum(struct photonList *photon_list);

void incrementNullPhotonNum(struct photonList *photon_list);

void verifyPhotonNum(struct photonList *photon_list);
