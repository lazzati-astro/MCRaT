
void readPlutoChombo( char pluto_file[200], double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);

void modifyPlutoName(char file[200], char prefix[200], int frame);
