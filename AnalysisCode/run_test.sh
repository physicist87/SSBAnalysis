#!/bin/bash

inputlists=("Data_SingleMuon_Run2016Bv2_1")
# inputlists=("DYJetsToLL_M_10To50_1")
runPeriod="UL2016PostVFP"
StudyName="Testv1"
Channels="MuMu"
Sample="TTbar_Signal"
echo $runPeriod
inputlists=("TTbar_Signal_1")
configdir=""
confch=""
echo "Good"

if [ "$Channels" = "MuMu" ]; then
    confch="dimuon.config"
elif [ "$Channels" = "ElEl" ]; then
    confch="dielec.config"
elif [ "$Channels" = "MuEl" ]; then
    confch="muelec.config"
else
    echo "Unknown Channels: $Channels"
    exit 1
fi

configpath="ULSummer20/${runPeriod}/"

for i in "${inputlists[@]}"; do
   mkdir -p output/${StudyName}/${runPeriod}/${Channels}/${Sample}
   # command line
   ./ssb_analysis ${runPeriod}/${Sample}/${i}.list ${StudyName}/${runPeriod}/${Channels}/$Sample/${i}.root ${configpath}${confch}
   echo ./ssb_analysis ${runPeriod}/${i}.list ${StudyName}/${runPeriod}/${Channels}/${i}.root ${configpath}${confch}
done

