
typedef struct box2d {
   int  lo_i;
   int  lo_j;
   int  hi_i;
    int hi_j;
} box2d;

void readPlutoChombo( char pluto_file[STR_BUFFER], double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);

void modifyPlutoName(char file[STR_BUFFER], char prefix[STR_BUFFER], int frame);

