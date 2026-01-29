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

struct photon createPhoton();

void initializePhoton(struct photon *ph);

double samplePhotonPhi(gsl_rng * rand, FILE *fPtr);

double samplePhotonTheta(double *velocity, gsl_rng * rand, FILE *fPtr);

void saveUserDefinePhoton(struct photon *ph_orig, struct photon *ph_user, struct hydro_dataframe *hydro_data, int hydro_index, gsl_rng * rand, FILE *fPtr);
