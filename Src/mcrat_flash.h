//
//  mcrat_flash.h
//  
//
//  Created by Tyler Parsotan on 7/23/21.
//

void modifyFlashName(char flash_file[STR_BUFFER], char prefix[STR_BUFFER], int frame);

void readAndDecimate(char flash_file[STR_BUFFER], struct hydro_dataframe *hydro_data, double r_inj, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr);
