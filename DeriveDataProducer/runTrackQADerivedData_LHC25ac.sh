#!/bin/bash

# env | grep prox
# unset http_proxy https_proxy
# env | grep prox

date
starttime=$(date)

#01- token check
alien-token-info &> tokenStatus.log
if [[ $(tail -1 tokenStatus.log) == "File >>><<< not found" ]];then  
   echo "ERROR :: Token Not Found, Create the token" 
   echo "$(date) ERROR :: Token Not Found, Create the token" >> tokenError.log ; for i in $(seq 1 10); do echo -en "\a"; sleep 0.25; done
   exit 1  #stop execution if token not found
elif [[ $(tail -1 tokenStatus.log) == "alien-token-info: command not found" ]]; then
   echo "ERROR :: alien-token-info command not found"
   echo "do ==> alienv enter O2Physics/latest ninja/latest" ; for i in $(seq 1 10); do echo -en "\a"; sleep 0.25; done
   exit 1  #stop execution if environment is not loaded
else
   echo "Token Found, Processing further ..."
fi

#02- Build The Task 
echo "BUILDING EXECUTABLE"
currentDir=$(pwd)
cd ~/alice/sw/BUILD/O2Physics-latest/O2Physics/
time ninja install &> output-01.log;
if [[ $(tail -1 output-01.log) == *Up-to-date* ]] ; then
   echo "BUILD SUCESSFULL"
else 
   if tail -1 output-01.log | grep -q "ninja: command not found"; then
   echo "ninja: command not found :: do ==> alienv enter O2Physics/latest ninja/latest"
   for i in $(seq 1 10); do echo -en "\a"; sleep 0.25; done
   exit 1 
   else
   echo "BUILD FAILED" ; code output-01.log ; for i in $(seq 1 10); do echo -en "\a"; sleep 0.25; done
   exit 1 #stop execution if build failed
   fi
fi

cd $currentDir
date
echo PROCESSING

rm -rf AnalysisResults.root
rm -rf AO2D_Derived.root

OutFile2="trkQADerive_output-02.log"
OutFile3="trkQADerive_output-03.log"

ROOT_File="Skimmed-LHC25ac_apass1_0.root"
 INPUT_JSON="LHC25ac_trackQADeriveInput_2025.05.23.json"
WRITER_JSON="LHC25ac_trackQADeriveOutput_2025.05.23.json"

export OPTIONS=" --configuration json://$INPUT_JSON  --resources-monitoring 2  --aod-memory-rate-limit 1000000000  --shm-segment-size 7500000000"


time o2-analysis-ft0-corrected-table   ${OPTIONS} |\
o2-analysis-hf-candidate-selector-d0 ${OPTIONS} |\
o2-analysis-track-to-collision-associator ${OPTIONS} |\
o2-analysis-hf-track-index-skim-creator ${OPTIONS} |\
o2-analysis-hf-pid-creator ${OPTIONS} |\
o2-analysis-hf-candidate-creator-2prong ${OPTIONS} |\
o2-analysis-pid-tof-beta          ${OPTIONS} |\
o2-analysis-trackselection        ${OPTIONS} |\
o2-analysis-pid-tof-full          ${OPTIONS} |\
o2-analysis-lf-lambdakzerobuilder ${OPTIONS} |\
o2-analysis-pid-tpc-base          ${OPTIONS} |\
o2-analysis-pid-tpc               ${OPTIONS} |\
o2-analysis-multiplicity-table    ${OPTIONS} |\
o2-analysis-event-selection       ${OPTIONS} |\
o2-analysis-pid-tof-base          ${OPTIONS} |\
o2-analysis-timestamp             ${OPTIONS} |\
o2-analysis-track-propagation     ${OPTIONS} |\
o2-analysis-occ-table-producer    ${OPTIONS} |\
o2-analysis-track-qa-derivedata   ${OPTIONS} --aod-file $ROOT_File --aod-writer-json $WRITER_JSON \
&> $OutFile2

echo $starttime >> $OutFile2
date >> $OutFile2
code $OutFile2

grep -e "IST" -e "\\[ERROR\\]" -e "\\[FATAL\\]" -e "segmentation" -e "Segmentation" -e "SEGMENTATION" -e "command not found" -e "Error:" -e "Error in " -e "\\[WARN\\]" -e "DEBUG" $OutFile2  >> $OutFile3

echo "Printing the Errors and Warnings " &> $OutFile3 
grep -e "IST" -e "\\[ERROR\\]" -e "\\[FATAL\\]" -e "segmentation" -e "Segmentation" -e "SEGMENTATION" -e "command not found" -e "Error:" -e "Error in " -e "\\[WARN\\]" -e "DEBUG" $OutFile2  >> $OutFile3
code $OutFile3
date

for i in $(seq 1 10); do echo -en "\a" ; sleep 0.25; done

echo "PROCESSING STATUS::COMPLETED"
