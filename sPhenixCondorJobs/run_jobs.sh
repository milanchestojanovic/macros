#! /bin/bash
 export USER="$(id -u -n)"
 export LOGNAME=${USER}
 export HOME=/sphenix/u/${LOGNAME}

 source /opt/sphenix/core/bin/sphenix_setup.sh -n

 # print the environment - needed for debugging
 printenv
 # this is how you run your Fun4All_G4_sPHENIX.C macro in batch: 


 # Arguments passed by Condor
 NUM_EVENTS=$1
 PROCESS_ID=$2
 INPUT_LIST_FILE=$3

 # Get the input file name for this process from the input list file
 INPUT_FILE=$(sed -n "${PROCESS_ID}p" ${INPUT_LIST_FILE})

 # Define the output file name
 OUTPUT_FILE=" /gpfs/mnt/gpfs02/sphenix/user/mstojano/output_${PROCESS_ID}.root"

 echo "Processing job "$2

 cd baseDir=/sphenix/u/mstojano/wsuwork/tutorials/CaloDataAnaRun24pp/macro

 # Run the ROOT macro with the specified arguments
 root.exe -l -b -q "Fun4All_CaloTreeGen.C(${NUM_EVENTS}, \"${INPUT_FILE}\", \"${OUTPUT_FILE}\")"

 echo all done
