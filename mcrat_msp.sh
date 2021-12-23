#!/bin/sh

#  mcrat_msp.sh
#  stands for mcrat (m)ost (s)cattering (p)rocesses
#
#  Shell script to determine which MCRaT processes have a very large number of scatterings and
#  allows the user to delete those porocesses chkpt and proc files to force a new restart with
#  less injected photons. Also allows user to refine mc.par file to set a new
#  number of injected photons for these same processs that will be restarted the next time the user
#  starts mcrat with the 'c' flag to continue a simulation.
#
#  Just calling mcrat_msp with an angle directory (e.g /path/to/hydro/and/MCRaT/0.0-3.0 ) and a number N_scatt
#  will find and list all the processes in the directory that have under gone > N_scatt scattering events.
#
#  Calling mcrat_msp with the -delete flag, an angle directory, and N_scatt will delete all of the chkpt and
#  proc files of the proceses in the angle directory that have undergone > N_scatt scattering events.
#
#  In order to modify the mc.par file to inject less photons call mcrat_msp with the -modifypar flag followed by
#  2 arguments: N_lower and N_higher for the lower and higher limits of the number of photons for mcrat to inject
#  initially. This assumes that the file structure is organized in the same manner outlined in the MCRaT manual.
#
#  All of the steps can be completed at once, however it is recommended that the user first list all of the files
#  that may be deleted by the script and verify them before running mcrat_msp to delete the files and modify the mc.par file.
#
#  Example calls are:
#  ./mcrat_msp.sh -modifypar 50 100 ~/path/to/mc_dir/0.0-3.0/ 1230000 -delete
#  ./mcrat_msp.sh ~/path/to/mc_dir/0.0-3.0/ 1230000
#
#  Created by Tyler Parsotan on 11/20/18.
#

RED='\033[0;31m'
NC='\033[0m' # No Color

DELETE_FLAG=0
MODIFY_FLAG=0

#read in the flags and directory given to script
while [ "$1" != "" ];
#echo "$1"
    do case $1 in
        -delete )
            DELETE_FLAG=1
        ;;
        -modifypar )
            shift; NMIN=$1 NMAX=$2 MODIFY_FLAG=1
        ;;
    esac;

    if [ -d "$1" ]
    then
       DIRECTORY=$1
        shift;
       N_SCATT_USER=$1
    fi

    shift;
done

#echo $DELETE_FLAG
#echo $MODIFY_FLAG
#echo $NMIN
#echo $NMAX
printf "In the directory %s:\n" "$DIRECTORY"
#echo $N_SCATT_USER
#  get files that aren't completed in directory: INCOMP_FILES= $(grep -L "completed" ~/Downloads/TEST_LOGS/mc_output_*)
INCOMP_FILES=$(grep -L "completed" "${DIRECTORY}"/mc_output_*)

#  go through files to se which has a number of scatterings larger than the specified N_scatt, if process hasn't even completed first frame, just look at the last scattering number
for i in $INCOMP_FILES;
do
    #  list the process numbers and number of scatterings

    FILE_NAME=${i##*/} #get the mc_output_NUM filename
    NO_END=${FILE_NAME%.*}
    PROC_NUM=${NO_END##*_} #get third value of filenme parsed by underscore
    NUM_FRAMES=$(grep "Working on frame" $i |wc -l)

    if (( "$NUM_FRAMES" > 1 )); #if the number of frames completed is >0
    then
        LINE=$(grep "Working on frame" $i | tail -n 2| head -n 1)

        NUM_SCATT_LINE=$(grep "The number of scatterings in this frame is:" $i |tail -n 1) # what if the process didnt even finish the first frame? then go to else
        #get the number of scatterings
        NUM_SCATT=${NUM_SCATT_LINE##*:}

        printf "Process %s recently completed frame: %s\n" "$PROC_NUM" "${LINE##* }"
        echo  $NUM_SCATT_LINE
    else
        LINE=$(grep "Working on frame" $i | tail -n 1)

        NUM_SCATT_LINE=$(grep "Scattering Number" $i |tail -n 1) #if the process didnt even finish the forst frame probably has alot alot of scatterings and should still have its files deleted in order to be restarted
        #get the number of scatterings
        NUM_SCATT=${NUM_SCATT_LINE##*:}
        printf "Process %s is working on the first frame number: %s\n" "$PROC_NUM" "${LINE##* }"
        printf "The number of scatterings completed in this frame is: %s\n" "$NUM_SCATT"

    fi



    if (("$NUM_SCATT" >= "$N_SCATT_USER"));
    then
        printf "${RED}Process %s has undergone more than %s scatterings${NC}\n" "$PROC_NUM" "$N_SCATT_USER"
        #  if the -delete flag is set, delete the chkpt and proc files
        # find ${directory} -name 'mc_chkpt_(number).dat' -delete #or 'mc_proc_(number).h5'
        if (("$DELETE_FLAG" == 1));
        then
            printf "${RED}Now deleting the checkpoint and process files${NC}\n"
            find "${DIRECTORY}" -name "mc_chkpt_${PROC_NUM}.dat" -delete
            find "${DIRECTORY}" -name "mc_proc_${PROC_NUM}.h5" -delete
        fi
    fi
    printf "\n"
#exit
done

#  if the -modifypar file is set modify the N_inj parameters in the mc.par file
if (("$MODIFY_FLAG" == 1));
then
    TAB=$'\t'
    #the min and max number of photons to inject is the 16th and 17th lines
    sed "21s/.*/$NMIN${TAB}${TAB}#Min number of photons/" "${DIRECTORY%/*/*}"/mc.par > "${DIRECTORY%/*/*}"/new_file.txt #replace nmin
    sed  -i '' -e "22s/.*/$NMAX${TAB}${TAB}#Max number of photons/" "${DIRECTORY%/*/*}"/new_file.txt #replace nmin

    mv "${DIRECTORY%/*/*}"/mc.par "${DIRECTORY%/*/*}"/old_mc.par
    mv "${DIRECTORY%/*/*}"/new_file.txt "${DIRECTORY%/*/*}"/mc.par #keep copy of old mc.par and move new_file (modified mc.par) to become mc.par
    printf "${RED}Make sure that the mc.par file is set to continue the simulation with 'c' and not 'i'${NC}\n"
fi

