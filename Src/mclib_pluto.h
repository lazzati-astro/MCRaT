
typedef struct box2d {
   int  lo_i;
   int  lo_j;
#if DIMENSIONS == THREE
   int  lo_k;
#endif
   int  hi_i;
   int hi_j;
#if DIMENSIONS == THREE
    int  hi_k;
#endif
} box2d;

typedef struct box3d {
   int  lo_i;
   int  lo_j;
   int  lo_k;
   int  hi_i;
   int hi_j;
    int  hi_k;
} box3d;

void readPlutoChombo( char pluto_file[STR_BUFFER], struct hydro_dataframe *hydro_data, double r_inj, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);

void modifyPlutoName(char file[STR_BUFFER], char prefix[STR_BUFFER], int frame);

void readGridFile(char gridfile[STR_BUFFER], double **grid_x1, double **grid_x2, double **grid_x3, double **grid_dx1, double **grid_dx2, double **grid_dx3, int *array_size, FILE *fPtr);

void readDblOutFile(char dblfile[STR_BUFFER], int *num_var, char ***var_strings, FILE *fPtr);

void readPluto(char pluto_file[STR_BUFFER], struct hydro_dataframe *hydro_data, double r_inj, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);
