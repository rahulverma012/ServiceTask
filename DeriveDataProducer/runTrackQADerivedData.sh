#!/bin/bash
# env | grep prox
# unset http_proxy https_proxy
# env | grep prox

date
starttime=$(date)
rm -rf AnalysisResults.root
rm -rf AO2D_Derived.root

OutFile2="Derive_output-02.log"
OutFile3="Derive_output-03.log"

currentDir=$(pwd)
cd ~/alice/sw/BUILD/O2Physics-latest/O2Physics
ninja install &> output-01.log ; code output-01.log 

cd $currentDir

ROOT_File="Skimmed-Run_55990n_AO2D.root"
# ROOT_File="Run_55990n_AO2D.root"
export OPTIONS=" --configuration json://myConfigDeriveData.json  --resources-monitoring 2  --aod-memory-rate-limit 1000000000  --shm-segment-size 7500000000"

time o2-analysis-ft0-corrected-table   ${OPTIONS} |\
o2-analysis-pid-tof-full          ${OPTIONS} |\
o2-analysis-pid-tpc-base          ${OPTIONS} |\
o2-analysis-pid-tpc               ${OPTIONS} |\
o2-analysis-occ-table-producer    ${OPTIONS} |\
o2-analysis-pid-tof-base          ${OPTIONS} |\
o2-analysis-pid-tof-beta          ${OPTIONS} |\
o2-analysis-event-selection       ${OPTIONS} |\
o2-analysis-multiplicity-table    ${OPTIONS} |\
o2-analysis-lf-lambdakzerobuilder ${OPTIONS} |\
o2-analysis-timestamp             ${OPTIONS} |\
o2-analysis-track-propagation     ${OPTIONS} |\
o2-analysis-track-qa-derivedata   ${OPTIONS} --aod-file $ROOT_File --aod-writer-json OutputDirectorDeriveData.json \
&> $OutFile2


# --aod-writer-ntfmerge 50
# --fairmq-ipc-prefix . \
# ->Buildoutput.log

# {
#     "OutputDirector": {
#         "debug_mode": true,
#         "resfile": "AO2D",
#         "OutputDescriptors": [
#             {
#                 "table": "AOD/COLLISION/0"
#             },
#             {
#                 "table": "AOD/MYTABLE/0"
#             }
#         ],
#         "ntfmerge": 1
#        original DF= 2336960459657024
#                              1000000
#                     3000000000000000
#                              1000000 =>     1M             1,000,000 - Gives two folders - 6993, 2334
#                             10000000 =>    10M            10,000,000 - Same as above
#                            100000000 =>   100M           100,000,000 - This Worked,  9327 = 6993+2334
#                           1000000000 =>     1B         1,000,000,000 - Working Properly ==> range of int is  upto 2,147,483,647 (maybe it is int)
#                          10000000000 =>    10B        10,000,000,000 - Working but folder name is different, DF_2336959213360128
#                         100000000000 =>   100B       100,000,000,000 - Working but folder name is different, DF_2336960199036928
#                        1000000000000 =>     1T     1,000,000,000,000 - It Failed
#                       10000000000000 =>    10T    10,000,000,000,000 - Worked but folder name is different   DF_2336959420850176
#                      100000000000000 =>   100T   100,000,000,000,000 - Worked but again different DF Name    DF_2336960221626368
#                     1000000000000000 =>  1000T 1,000,000,000,000,000 - It Failed
#     }
# }

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
