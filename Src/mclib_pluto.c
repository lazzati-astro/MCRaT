//
//  mclib_pluto.c
//  to read in pluto AMR chombo files
//
//  Created by Tyler Parsotan on 11/26/19.
//

#include "mcrat.h"

#define PLUTO_VAR_LENGTH 5

void readPlutoChombo( char pluto_file[STR_BUFFER], struct hydro_dataframe *hydro_data, double r_inj, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr)
{
    hid_t  file, dset, space, group, attr, filetype, memtype;
    herr_t status;
    hsize_t dims[1]={0}; //hold number of processes in each level
    size_t  sdim=0;
    int i=0, j=0, k=0, l=0, m=0, n=0, r_count=0, num_dims=0, num_levels=0, num_vars=4, logr=0;
    int nbx=0, nby=0, nbz=0, total_size=0, total_box_size=0, offset=0, elem_factor=0;
    #if DIMENSIONS == THREE
        box2d prob_domain[1]={0,0,0,0,0,0};
    #else
        box2d prob_domain[1]={0,0,0,0};
    #endif
    char level[200]="", component_num[200]="";
    double ph_rmin=0, ph_rmax=0, ph_thetamin=0, ph_thetamax=0, r_grid_innercorner=0, r_grid_outercorner=0, theta_grid_innercorner=0, theta_grid_outercorner=0;
    int *level_dims=NULL, *box_offset=NULL;
    double *x1_array=NULL, *x2_array=NULL, *dx1_array=NULL, *dx2_array=NULL, *x3_array=NULL, *dx3_array=NULL;
    double *dombeg1=NULL, *dombeg2=NULL, *dombeg3=NULL, *dx=NULL, *g_x2stretch=NULL, *g_x3stretch=NULL;
    double *all_data=NULL, *x1_buffer=NULL, *x2_buffer=NULL, *x3_buffer=NULL, *dx1_buffer=NULL, *dx2_buffer=NULL, *dx3_buffer=NULL;
    double *dens_buffer=NULL, *pres_buffer=NULL, *vel_x2_buffer=NULL, *vel_x1_buffer=NULL, *vel_x3_buffer=NULL;
    double *B_x2_buffer=NULL, *B_x1_buffer=NULL, *B_x3_buffer=NULL;
    char **var_strings=NULL;
    
    if (ph_inj_switch==0)
    {
        ph_rmin=min_r;
        ph_rmax=max_r;
        ph_thetamin=min_theta-2*0.017453292519943295; //min_theta - 2*Pi/180 (2 degrees)
        ph_thetamax=max_theta+2*0.017453292519943295; //max_theta + 2*Pi/180 (2 degrees)
    }

    
    //define dataset for boxes of each level
    hid_t box_dtype = H5Tcreate (H5T_COMPOUND, sizeof(box2d));
    H5Tinsert(box_dtype, "lo_i", HOFFSET(box2d, lo_i), H5T_NATIVE_INT);
    H5Tinsert(box_dtype, "lo_j", HOFFSET(box2d, lo_j), H5T_NATIVE_INT);
    #if DIMENSIONS == THREE
        H5Tinsert(box_dtype, "lo_k", HOFFSET(box2d, lo_k), H5T_NATIVE_INT);
    #endif
    H5Tinsert(box_dtype, "hi_i", HOFFSET(box2d, hi_i), H5T_NATIVE_INT);
    H5Tinsert(box_dtype, "hi_j", HOFFSET(box2d, hi_j), H5T_NATIVE_INT);
    #if DIMENSIONS == THREE
        H5Tinsert(box_dtype, "hi_k", HOFFSET(box2d, hi_k), H5T_NATIVE_INT);
    #endif
    
    //open the pluto file
    file = H5Fopen (pluto_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    fprintf(fPtr, ">> MCRaT: Reading positional, density, pressure, and velocity information...\n");
    fflush(fPtr);

    //1. read in the number of dims
    group = H5Gopen (file, "/Chombo_global", H5P_DEFAULT);
    
    attr = H5Aopen (group, "SpaceDim", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_dims);
    
    status = H5Aclose (attr);
    status = H5Gclose (group);
    
    //printf("readPlutoChombo spacedim: %d\n", num_dims);

    
    //2. get the number of levels
    attr = H5Aopen (file, "num_levels", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_levels);
    
    status = H5Aclose (attr);
    //printf("readPlutoChombo num_levels: %d\n", num_levels);
    
    dombeg1=malloc(num_levels*sizeof(double));
    dombeg2=malloc(num_levels*sizeof(double));
    dombeg3=malloc(num_levels*sizeof(double));
    dx=malloc(num_levels*sizeof(double));
    g_x2stretch=malloc(num_levels*sizeof(double));
    g_x3stretch=malloc(num_levels*sizeof(double));
    level_dims=malloc(num_levels*sizeof(int));
    
    //3. get number of variables to read in (should be 4 in non-MHD case) and read in variable type and order
    attr = H5Aopen (file, "num_components", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_vars);
    status = H5Aclose (attr);
    
    var_strings = malloc(num_vars * sizeof(char*));
    for (i=0;i<num_vars;i++)
    {
        //create name of component in hdf5 and malloc the memory to read variable string name
        snprintf(component_num,sizeof(component_num), "component_%d",i);
        
        attr = H5Aopen (file, component_num, H5P_DEFAULT);
        
        //Get the datatype and its size and allocate memory
        filetype = H5Aget_type (attr);
        sdim = H5Tget_size (filetype);
        sdim++;                         /* Make room for null terminator */
        var_strings[i] = malloc(sdim * sizeof(char));


        //Create the memory datatype.
        memtype = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (memtype, sdim);

        //Read the data.
        status = H5Aread (attr, memtype, var_strings[i]);
        
        status = H5Aclose (attr);
        status = H5Tclose (filetype);
        status = H5Tclose (memtype);
    }
    
    //printf("readPlutoChombo num_vars: %d\n", num_vars);
    
    //get the total number of values that I need to allocate memory for
    for (i=0;i<num_levels;i++)
    {
        snprintf(level, sizeof(level), "level_%d", i);
        //printf("Opening level %d Boxes\n", i);
        
        group = H5Gopen(file, level, H5P_DEFAULT);
        
        dset= H5Dopen(group, "data:datatype=0", H5P_DEFAULT);
        
        //get dimensions of array and save it
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
        
        total_size+=dims[0];
        *(level_dims+i)=dims[0];
        
        status = H5Sclose (space);
        H5Dclose(dset);
        H5Gclose(group);
    }
    //printf("The total number of elements is %d\n", total_size);
    
    //now allocate space to save data
    all_data=malloc(total_size*sizeof (double));
    
    //and allocate arrays to hold r and theta values to calculate them on the fly
    x1_buffer=malloc((total_size/num_vars)*sizeof (double));
    x2_buffer=malloc((total_size/num_vars)*sizeof (double));
    dx1_buffer=malloc((total_size/num_vars)*sizeof (double));
    dx2_buffer=malloc((total_size/num_vars)*sizeof (double));

    vel_x1_buffer= malloc ((total_size/num_vars) * sizeof (double));
    vel_x2_buffer=malloc ((total_size/num_vars) * sizeof (double));
    dens_buffer= malloc ((total_size/num_vars) * sizeof (double));
    pres_buffer=malloc ((total_size/num_vars) * sizeof (double));
    
//NEED TO IMPLEMENT THIS BELOW
    #if B_FIELD_CALC == SIMULATION
        B_x1_buffer= malloc ((total_size/num_vars) * sizeof (double));
        B_x2_buffer=malloc ((total_size/num_vars) * sizeof (double));
    #endif

    
    #if DIMENSIONS == THREE
        x3_buffer=malloc((total_size/num_vars)*sizeof (double));
        dx3_buffer=malloc((total_size/num_vars)*sizeof (double));
    #endif
    
    #if DIMENSIONS == THREE || DIMENSIONS == TWO_POINT_FIVE
        vel_x3_buffer=malloc ((total_size/num_vars) * sizeof (double));
        #if B_FIELD_CALC==SIMULATION
            B_x3_buffer= malloc ((total_size/num_vars) * sizeof (double));
        #endif
    #endif
    
    offset=0;
    //read in the data
    for (i=0;i<num_levels;i++)
    {
        snprintf(level, sizeof(level), "level_%d", i);
        //printf("Opening level %d Boxes\n", i);
        
        group = H5Gopen(file, level, H5P_DEFAULT);
        
        //read in the data
        dset= H5Dopen(group, "data:datatype=0", H5P_DEFAULT);
        status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,(all_data+offset));
        H5Dclose(dset);
        //printf("First few data %e %e %e Last few data %e %e\n", *(all_data+offset), *(all_data+offset+1), *(all_data+offset+2), *(all_data+offset+(*(level_dims+i))-2), *(all_data+offset+(*(level_dims+i))-1));
        
        //read in the box offsets in all_data
        dset= H5Dopen(group, "data:offsets=0", H5P_DEFAULT);
        
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
        box_offset=malloc(dims[0]*sizeof(int));
        status = H5Sclose (space);
        
        status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_offset);
        H5Dclose(dset);

        //open various properties about that refinement level
        attr = H5Aopen (group, "prob_domain", H5P_DEFAULT);
        status = H5Aread (attr, box_dtype, &prob_domain);
        status = H5Aclose (attr);
        //printf("Prob_domain %d %d %d %d\n", prob_domain->lo_i, prob_domain->lo_j, prob_domain->hi_i, prob_domain->hi_j);
               
        attr = H5Aopen (group, "dx", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dx+i));
        status = H5Aclose (attr);
        //printf("dx %e\n", *(dx+i));
        
        attr = H5Aopen (group, "logr", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_INT, &logr);
        status = H5Aclose (attr);
        //printf("logr %d\n", logr);
        
        attr = H5Aopen (group, "domBeg1", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dombeg1+i));
        status = H5Aclose (attr);
        //printf("dombeg1 %e\n", *(dombeg1+i));
        
        //set default just in case
        *(dombeg2+i)=0;
        *(dombeg3+i)=0;
        
        attr = H5Aopen (group, "g_x2stretch", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_DOUBLE, (g_x2stretch+i));
        status = H5Aclose (attr);
        //printf("g_x2stretch %e\n", *(g_x2stretch+i));
        
        attr = H5Aopen (group, "domBeg2", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dombeg2+i));
        status = H5Aclose (attr);
        //printf("dombeg2 %e\n", *(dombeg2+i));

        #if DIMENSIONS == THREE
            attr = H5Aopen (group, "g_x3stretch", H5P_DEFAULT);
            status = H5Aread (attr, H5T_NATIVE_DOUBLE, (g_x3stretch+i));
            status = H5Aclose (attr);
            
            attr = H5Aopen (group, "domBeg3", H5P_DEFAULT);
            status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dombeg3+i));
            status = H5Aclose (attr);
        #endif
               
        //read in the boxes
        dset= H5Dopen(group, "boxes", H5P_DEFAULT);
        
        //get dimensions of array and save it
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
                    
        box2d box_data[dims[0]];
        
        status = H5Dread (dset, box_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,box_data);
        status = H5Sclose (space);
        H5Dclose(dset);
        
        x1_array=malloc( (prob_domain->hi_i - prob_domain->lo_i +1) * sizeof (double ));
        x2_array=malloc( (prob_domain->hi_j - prob_domain->lo_j +1) * sizeof (double ));
        dx1_array=malloc( (prob_domain->hi_i - prob_domain->lo_i +1) * sizeof (double ));
        dx2_array=malloc( (prob_domain->hi_j - prob_domain->lo_j +1) * sizeof (double ));
        #if DIMENSIONS == THREE
            x3_array=malloc( (prob_domain->hi_k - prob_domain->lo_k +1) * sizeof (double ));
            dx3_array=malloc( (prob_domain->hi_k - prob_domain->lo_k +1) * sizeof (double ));
        #endif

        //create the arrays that hold the refinement level radii and angles
        for (j=0;j<(prob_domain->hi_i - prob_domain->lo_i +1);j++)
        {
            if (logr==0)
            {
                *(x1_array+j)=(*(dombeg1+i)) + (*(dx+i)) * (prob_domain->lo_i + j + 0.5);
                *(dx1_array+j)=(*(dx+i));
            }
            else
            {
                *(x1_array+j)=(*(dombeg1+i)) * 0.5 * (exp((*(dx+i)) * (prob_domain->lo_i + j + 1)) + exp((*(dx+i)) * (prob_domain->lo_i + j))   );
                *(dx1_array+j)=(*(dombeg1+i)) * (exp((*(dx+i)) * (prob_domain->lo_i + j + 1)) - exp((*(dx+i)) * (prob_domain->lo_i + j))   );
            }
            //if (i==2) printf("x1_array: %0.8e dr: %0.8e\n", *(x1_array+j), *(dx1_array+j));
        }
        
        for (j=0;j<(prob_domain->hi_j - prob_domain->lo_j +1);j++)
        {
            *(x2_array+j)=(*(dombeg2+i)) + (*(dx+i)) * (*(g_x2stretch+i)) * (prob_domain->lo_j + j + 0.5);
            *(dx2_array+j)=(*(dx+i))*(*(g_x2stretch+i));
        }
        
        #if DIMENSIONS == THREE
            for (j=0;j<(prob_domain->hi_k - prob_domain->lo_k +1);j++)
            {
                *(x3_array+j)=(*(dombeg3+i)) + (*(dx+i)) * (*(g_x3stretch+i)) * (prob_domain->lo_k + j + 0.5);
                *(dx3_array+j)=(*(dx+i))*(*(g_x3stretch+i));
            }
        #endif
        
        //go through the boxes to create the buffer arrays
        total_box_size=0;
        for (j=0; j<dims[0]; j++)
        {
            //printf("i %d %d %d %d %d \n", j, box_data[j].lo_i, box_data[j].lo_j, box_data[j].hi_i, box_data[j].hi_j);
            nbx=box_data[j].hi_i-box_data[j].lo_i+1;
            nby=box_data[j].hi_j-box_data[j].lo_j+1;
            nbz=1;
            
            #if DIMENSIONS == THREE
                nbz=box_data[j].hi_k-box_data[j].lo_k+1;
            #endif
            
            //loop over each variable values of box
            for (k=0;k<num_vars;k++)
            {
                //loop over the x1, x2, x3 variables
                for (l=0; l<nbz; l++)
                {
                    //loop over the theta
                    for (m=0 ;m<nby ;m++)
                    {
                        //loop over radii
                        for (n=0; n< nbx;n++)
                        {
                        
                            *(x1_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx*nby + m*nbx + n )= (*(x1_array+box_data[j].lo_i+n));
                            *(x2_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n )=(*(x2_array+box_data[j].lo_j+m));
                            *(dx1_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(dx1_array+box_data[j].lo_i+n));
                            *(dx2_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n )=(*(dx2_array+box_data[j].lo_j+m));
                            #if DIMENSIONS == THREE
    //NEED TO FIGURE OUT INDEXES HERE and elsewhere in this loop for the 3D case
                                *(x3_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n )=(*(x3_array+box_data[j].lo_k+l));
                                *(dx3_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n )=(*(dx3_array+box_data[j].lo_k+l));
                            #endif
                            
                            //multipy by hydro caling factors
                            *(x1_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx*nby + m*nbx + n ) *= HYDRO_L_SCALE;
                            *(dx1_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n) *= HYDRO_L_SCALE;
                            
                            #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                                *(x2_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx*nby + m*nbx + n  ) *= HYDRO_L_SCALE;
                                *(dx2_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n) *= HYDRO_L_SCALE;
                            #endif
                            
                            #if DIMENSIONS == THREE
                                #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                                    *(x3_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx*nby + m*nbx + n) *= HYDRO_L_SCALE;
                                    *(dx3_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n) *= HYDRO_L_SCALE;
                                #endif
                            #endif

                            
                            if (strcmp(var_strings[k], "rho") == 0)
                            {
                                *(dens_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n))*HYDRO_D_SCALE;
                            }
                            else if (strcmp(var_strings[k], "vx1") == 0)
                            {
                                *(vel_x1_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ));
                            }
                            else if (strcmp(var_strings[k], "vx2") == 0)
                            {
                                *(vel_x2_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ));
                            }
                            else if (strcmp(var_strings[k], "prs") == 0)
                            {
                                *(pres_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ))*HYDRO_P_SCALE;
                            }
                            #if B_FIELD_CALC == SIMULATION
                                else if (strcmp(var_strings[k], "bx1") == 0)
                                {
                                    *(B_x1_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n) = (*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ))*HYDRO_B_SCALE;
                                }
                                else if (strcmp(var_strings[k], "bx2") == 0)
                                {
                                    *(B_x2_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ))*HYDRO_B_SCALE;
                                }
                            #endif
                            #if DIMENSIONS == THREE || DIMENSIONS == TWO_POINT_FIVE
                                else if (strcmp(var_strings[k], "vx3") == 0)
                                {
                                    *(vel_x3_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ));
                                }
                                #if B_FIELD_CALC==SIMULATION
                                    else if (strcmp(var_strings[k], "bx3") == 0)
                                    {
                                        *(B_x3_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx*nby + m*nbx + n)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby*nbz + l*nbx*nby + m*nbx + n  ))*HYDRO_B_SCALE;
                                    }
                                #endif
                            #endif
                            
                        /*
                        switch (k) //instead compare string from component string to 'prs', 'vx1', etc
                        {
                            case 0:
                                //when k is 0 save the radii, caus eonly need to do once and also save the densities
                                *(x1_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx +m  )= (*(x1_array+box_data[j].lo_i+m));
                                *(x2_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx +m )=(*(x2_array+box_data[j].lo_j+l));
                                *(dx1_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx +m)=(*(dx1_array+box_data[j].lo_i+m));
                                *(dx2_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx +m )=(*(dx2_array+box_data[j].lo_j+l));
                                #if DIMENSIONS == THREE
//NEED TO FIGURE OUT INDEXES HERE and elsewhere in this loop for the 3D case
                                    *(x3_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx +m )=(*(x3_array+box_data[j].lo_j+l));
                                    *(dx3_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx +m )=(*(dx3_array+box_data[j].lo_j+l));
                                #endif
                                
                                //multipy by hydro caling factors
                                *(x1_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx +m  ) *= HYDRO_L_SCALE;
                                *(dx1_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx +m) *= HYDRO_L_SCALE;
                                
                                #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
                                    *(x2_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx +m  ) *= HYDRO_L_SCALE;
                                    *(dx2_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx +m) *= HYDRO_L_SCALE;
                                #endif
                                
                                #if DIMENSIONS == THREE
                                    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
                                        *(x3_buffer+offset/num_vars+(*(box_offset+j))/num_vars+ l*nbx +m  ) *= HYDRO_L_SCALE;
                                        *(dx3_buffer+(offset+(*(box_offset+j)))/num_vars + l*nbx +m) *= HYDRO_L_SCALE;
                                    #endif
                                #endif
                                
                                *(dens_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx +m)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby + l*nbx +m  ))*HYDRO_D_SCALE;
                                break;
                            case 1:
                                *(vel_x1_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx +m)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby + l*nbx +m  ));
                                break;
                            case 2:
                                *(vel_x2_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx +m)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby + l*nbx +m  ));
                                break;
                            case 3:
                                *(pres_buffer+ (offset+(*(box_offset+j)))/num_vars + l*nbx +m)=(*(all_data+offset+(*(box_offset+j))+ k*nbx*nby + l*nbx +m  ))*HYDRO_P_SCALE;
                                break;
                            default:
                                break;
                        }
                        */
                        }
                    }
                }
            }
            /*
              // this was for testing, the actual nested loop is above
              if (i==2 && j==dims[0]-1)
              {
                  
                  r_count=0;
                  for (k=0;k<1;k++)
                  {
                      //loop over the radii
                      for (l=0; l<nby; l++)
                      {
                          //loop over the angles
                          for (m=0 ;m<nbx ;m++)
                          {
                              printf("count: %d idx: %d all_data val: %0.8e r: %e theta %e\n", (*(box_offset+j))+r_count, (*(box_offset+j))+k*nbx*nby + l*nby +m, *(all_data+offset+(*(box_offset+j))+ k*nbx*nby + l*nbx +m  ), *(radii+box_data[j].lo_i+m), *(angles+box_data[j].lo_j+l));
                              r_count++;
                          }
                          printf("\n");
                      }
                      printf("\n");
                  }
                  
                  for (k=0;k<64;k++)
                  {
                      printf("r: %e, theta %e dens data: %e\n",*(x1_buffer+ (offset+(*(box_offset+j)))/num_vars+k), *(x2_buffer+ (offset+(*(box_offset+j)))/num_vars +k), *(dens_buffer+ (offset+(*(box_offset+j)))/num_vars +k) );
                  }
                  printf("\n");
                   
              }
        */
            
             
        }
        
        offset+=(*(level_dims+i));
        
        H5Gclose(group);
        free(x1_array); free(x2_array); free(x3_array); free(dx1_array); free(dx2_array); free(dx3_array); free(box_offset);
        
    }
    status = H5Fclose (file);
    free(dombeg1); free(dombeg2); free(dombeg3); free(dx); free(g_x2stretch); free(g_x3stretch); free(all_data);
    
    for (i=0;i<num_vars;i++)
    {
        free(var_strings[i]);
    }
    free(var_strings);
    
    //have all the data so need to go through and decide which values we will be keeping based on phtoon ranges we want to keep
    //fill in radius array and find in how many places r > injection radius
    #if CYCLOSYNCHROTRON_SWITCH == ON
        elem_factor=2;
    #else
        elem_factor=0;
    #endif
    r_count=0;
    while (r_count==0)
    {
        r_count=0;
        elem_factor++;
        for (i=0;i<(total_size/num_vars);i++)
        {
            if (ph_inj_switch==0)
            {
                #if DIMENSIONS == THREE
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i))-0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))-0.5*((*(dx2_buffer+i))), (*(x3_buffer+i))-0.5*(*(dx3_buffer+i)));
                        hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (*(x1_buffer+i))+0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))+0.5*((*(dx2_buffer+i))), (*(x3_buffer+i))+0.5*(*(dx3_buffer+i)));
                #else
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i))-0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))-0.5*((*(dx2_buffer+i))), 0);
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (*(x1_buffer+i))+0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))+0.5*((*(dx2_buffer+i))), 0);
                #endif
                    
                if (((ph_rmin - elem_factor*C_LIGHT/hydro_data->fps) <= r_grid_outercorner) && (r_grid_innercorner  <= (ph_rmax + elem_factor*C_LIGHT/hydro_data->fps) ) && (theta_grid_outercorner >= ph_thetamin) && (theta_grid_innercorner <= ph_thetamax) )
                {
                    r_count++;
                }
            
            }
            else
            {
                #if DIMENSIONS == THREE
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i)), (*(x2_buffer+i)), (*(x3_buffer+i)) );
                #else
                    hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i)), (*(x2_buffer+i)), 0);
                #endif

                if ( r_grid_innercorner > (0.95*r_inj) )
                {
                    r_count++;
                }
            }
        }
        //printf( "r_count: %d\n", r_count);
    }
    fprintf(fPtr, "Elem factor: %d Ph_rmin: %e rmax: %e Chosen FLASH min_r: %e max_r: %e min_theta: %e degrees max_theta: %e degrees\n", elem_factor, ph_rmin, ph_rmax, ph_rmin - (elem_factor*C_LIGHT/hydro_data->fps), ph_rmax + (elem_factor*C_LIGHT/hydro_data->fps), ph_thetamin*180/M_PI, ph_thetamax*180/M_PI);
    fflush(fPtr);

    
    //allocate memory to hold processed data
   (hydro_data->pres)=malloc (r_count * sizeof (double ));
   (hydro_data->v0)=malloc (r_count * sizeof (double ));
   (hydro_data->v1)=malloc (r_count * sizeof (double ));
   (hydro_data->dens)=malloc (r_count * sizeof (double ));
   (hydro_data->r0)=malloc (r_count * sizeof (double ));
   (hydro_data->r1)=malloc (r_count * sizeof (double ));
   (hydro_data->r)=malloc (r_count * sizeof (double ));
   (hydro_data->theta)=malloc (r_count * sizeof (double ));
   (hydro_data->gamma)=malloc (r_count * sizeof (double ));
   (hydro_data->dens_lab)=malloc (r_count * sizeof (double ));
   (hydro_data->r0_size)=malloc (r_count * sizeof (double ));
   (hydro_data->r1_size)=malloc (r_count * sizeof (double ));
   (hydro_data->temp)=malloc (r_count * sizeof (double ));
  
    #if B_FIELD_CALC == SIMULATION
       (hydro_data->B0)= malloc (r_count * sizeof (double));
       (hydro_data->B1)= malloc (r_count * sizeof (double));
    #endif


    #if DIMENSIONS == THREE
       (hydro_data->r2)=malloc(r_count*sizeof (double));
       (hydro_data->r2_size)=malloc(r_count*sizeof (double));
    #endif
                                               
    #if DIMENSIONS == THREE || DIMENSIONS == TWO_POINT_FIVE
       (hydro_data->v2)=malloc (r_count * sizeof (double));
        #if B_FIELD_CALC==SIMULATION
           (hydro_data->B2)= malloc (r_count * sizeof (double));
        #endif
    #endif

    
    //assign values based on r> 0.95*r_inj
    j=0;
    for (i=0;i<(total_size/num_vars);i++)
    {
        if (ph_inj_switch==0)
        {
            #if DIMENSIONS == THREE
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i))-0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))-0.5*((*(dx2_buffer+i))), (*(x3_buffer+i))-0.5*(*(dx3_buffer+i)));
                    hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (*(x1_buffer+i))+0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))+0.5*((*(dx2_buffer+i))), (*(x3_buffer+i))+0.5*(*(dx3_buffer+i)));
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i))-0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))-0.5*((*(dx2_buffer+i))), 0);
                hydroCoordinateToSpherical(&r_grid_outercorner, &theta_grid_outercorner, (*(x1_buffer+i))+0.5*(*(dx1_buffer+i)), (*(x2_buffer+i))+0.5*((*(dx2_buffer+i))), 0);
            #endif

            
            if (((ph_rmin - elem_factor*C_LIGHT/hydro_data->fps) <= r_grid_outercorner) && (r_grid_innercorner  <= (ph_rmax + elem_factor*C_LIGHT/hydro_data->fps) ) && (theta_grid_outercorner >= ph_thetamin) && (theta_grid_innercorner <= ph_thetamax))
            {
                ((hydro_data->pres))[j]=(*(pres_buffer+i));
                ((hydro_data->v0))[j]= (*(vel_x1_buffer+i));//*sin((*(x2_buffer+i))) + (*(vel_x2_buffer+i))*cos((*(x2_buffer+i))) ; //convert from spherical to cartesian
                ((hydro_data->v1))[j]= (*(vel_x2_buffer+i));//*cos((*(x2_buffer+i))) - (*(vel_x2_buffer+i))*sin((*(x2_buffer+i))); //z axis is actually y axis in 2D
            
                ((hydro_data->dens))[j]=(*(dens_buffer+i));
                ((hydro_data->r0))[j]=(*(x1_buffer+i)); //*sin((*(x2_buffer+i)));
                ((hydro_data->r1))[j]=(*(x2_buffer+i));//*cos((*(x2_buffer+i)));
                ((hydro_data->r))[j]=(*(x1_buffer+i));
                ((hydro_data->theta))[j]=(*(x2_buffer+i));//theta in radians in relation to jet axis,will be written over

                ((hydro_data->r0_size))[j]=(*(dx1_buffer+i));
                ((hydro_data->r1_size))[j]=(*(dx2_buffer+i));
                ((hydro_data->gamma))[j]=1/pow(1.0-( (*(vel_x1_buffer+i))*(*(vel_x1_buffer+i)) + (*(vel_x2_buffer+i))*(*(vel_x2_buffer+i)) ),0.5); //v is in units of c
                ((hydro_data->dens_lab))[j]= (*(dens_buffer+i)) /pow(1.0-( (*(vel_x1_buffer+i))*(*(vel_x1_buffer+i)) + (*(vel_x2_buffer+i))*(*(vel_x2_buffer+i)) ),0.5);
                ((hydro_data->temp))[j]=pow(3*(*(pres_buffer+i))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                
                #if B_FIELD_CALC == SIMULATION
                    (hydro_data->B0)[j]= (*(B_x1_buffer+i));
                    (hydro_data->B1)[j]= (*(B_x2_buffer+i));
                #endif

                #if DIMENSIONS == THREE
                    (hydro_data->r2)[j]=(*(x3_buffer+i));
                    (hydro_data->r2_size)[j]=(*(dx3_buffer+i));
                #endif
                                                           
                #if DIMENSIONS == THREE || DIMENSIONS == TWO_POINT_FIVE
                    (hydro_data->v2)[j]=(*(vel_x3_buffer+i));
                    #if B_FIELD_CALC==SIMULATION
                       (hydro_data->B2)[j]= (*(B_x3_buffer+i));
                    #endif
                #endif
                
                j++;
            }
        }
        else
        {
            #if DIMENSIONS == THREE
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i)), (*(x2_buffer+i)), (*(x3_buffer+i)) );
            #else
                hydroCoordinateToSpherical(&r_grid_innercorner, &theta_grid_innercorner, (*(x1_buffer+i)), (*(x2_buffer+i)), 0);
            #endif

            if ( r_grid_innercorner > (0.95*r_inj) )
            {
                ((hydro_data->pres))[j]=(*(pres_buffer+i));
                ((hydro_data->v0))[j]= (*(vel_x1_buffer+i));//*sin((*(x2_buffer+i))) + (*(vel_x2_buffer+i))*cos((*(x2_buffer+i))) ; //convert from spherical to cartesian
                ((hydro_data->v1))[j]= (*(vel_x2_buffer+i));//*cos((*(x2_buffer+i))) - (*(vel_x2_buffer+i))*sin((*(x2_buffer+i))); //z axis is actually y axis in 2D
            
                ((hydro_data->dens))[j]=(*(dens_buffer+i));
                ((hydro_data->r0))[j]=(*(x1_buffer+i)); //*sin((*(x2_buffer+i)));
                ((hydro_data->r1))[j]=(*(x2_buffer+i));//*cos((*(x2_buffer+i)));
                ((hydro_data->r))[j]=(*(x1_buffer+i));
                ((hydro_data->theta))[j]=(*(x2_buffer+i));//theta in radians in relation to jet axis, will be written over

                ((hydro_data->r0_size))[j]=(*(dx1_buffer+i));
                ((hydro_data->r1_size))[j]=(*(dx2_buffer+i));
                ((hydro_data->gamma))[j]=1/pow(1.0-( (*(vel_x1_buffer+i))*(*(vel_x1_buffer+i)) + (*(vel_x2_buffer+i))*(*(vel_x2_buffer+i)) ),0.5); //v is in units of c
                ((hydro_data->dens_lab))[j]= (*(dens_buffer+i)) /pow(1.0-( (*(vel_x1_buffer+i))*(*(vel_x1_buffer+i)) + (*(vel_x2_buffer+i))*(*(vel_x2_buffer+i)) ),0.5);
                ((hydro_data->temp))[j]=pow(3*(*(pres_buffer+i))*pow(C_LIGHT,2.0)/(A_RAD) ,1.0/4.0);
                    
                #if B_FIELD_CALC == SIMULATION
                    (hydro_data->B0)[j]= (*(B_x1_buffer+i));
                    (hydro_data->B1)[j]= (*(B_x2_buffer+i));
                #endif

                #if DIMENSIONS == THREE
                    (hydro_data->r2)[j]=(*(x3_buffer+i));
                    (hydro_data->r2_size)[j]=(*(dx3_buffer+i));
                #endif
                                                           
                #if DIMENSIONS == THREE || DIMENSIONS == TWO_POINT_FIVE
                    (hydro_data->v2)[j]=(*(vel_x3_buffer+i));
                    #if B_FIELD_CALC==SIMULATION
                       (hydro_data->B2)[j]= (*(B_x3_buffer+i));
                    #endif
                #endif

                j++;
            }
        }
    }
    
    //fprintf(fPtr, "number: %d\n", r_count);
    

   hydro_data->num_elements=r_count;
    
    
    free(x1_buffer); free(x2_buffer); free(x3_buffer); free(dx1_buffer); free(dx2_buffer); free(dx3_buffer); free(dens_buffer); free(pres_buffer);
    free(vel_x1_buffer); free(vel_x2_buffer); free(vel_x3_buffer); free(B_x1_buffer); free(B_x2_buffer); free(B_x3_buffer);

}

void modifyPlutoName(char file[STR_BUFFER], char prefix[STR_BUFFER], int frame)
{
    int lim1=0, lim2=0, lim3=0;
    char test[STR_BUFFER]="" ;
    
    #if DIMENSIONS == TWO || DIMENSIONS == TWO_POINT_FIVE
    {
        //2D case
        lim1=10;
        lim2=100;
        lim3=1000;
    }
    #else
    {
        //3d case
        lim1=100;
        lim2=1000;
        lim3=10000;
    }
    #endif
    
    if (frame<lim1)
    {
        snprintf(test,sizeof(test), "%s%.3d%d.hdf5",prefix,000,frame);
    }
    else if (frame<lim2)
    {
        snprintf(test,sizeof(test), "%s%.2d%d.hdf5",prefix,00,frame);
    }
    else if (frame<lim3)
    {
        snprintf(test,sizeof(test), "%s%d%d.hdf5",prefix,0,frame);
    }
    else
    {
        snprintf(test,sizeof(test), "%s%d.hdf5",prefix,frame);
    }
    strncpy(file, test, sizeof(test));//had to do this workaround for some weird reason
}

