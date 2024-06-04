#!/bin/bash

set inputlists = ("Data_DoubleMuon_Run2016B_193")
#set inputlists = ("Data_DoubleMuon_Run2016G_1")
#set inputlists = ("TTJets_Signal_GluoneMoveCRTune_erdON_1")
#set inputlists = ("TTJets_others_GluoneMoveCRTune_erdON_1")
#set inputlists = ("WJetsToLNu_1")
#set inputlists = ("TTJets_Signal_1")
set inputlists = ("Data_SingleMuon_Run2016Bv2_1")
#set inputlists = ("DYJetsToLL_M_10To50_1")
set runPeriod = "UL2016PostVFP"
set StudyName = "Testv1"
set Channels = "MuMu"
echo $runPeriod
set inputlists = ("TTbar_Signal_1")
set configdir = ""
set configch = ""
echo "Good"

if [ "$Channels" = "MuMu" ]; then
    confch="dimuon.cfg"
elif [ "$Channels" = "ElEl" ]; then
    confch="dielec.cfg"
elif [ "$Channels" = "MuEl" ]; then
    confch="muelec.cfg"
else
    echo "Unknown Channels: $Channels"
    exit 1
fi

set config = ULSummer20/${runPeriod}/

foreach i ( $inputlists )
   mkdir -p output/${StudyName}/${runPeriod}/${Channels}
   #./ssb_analysis ${runPeriod}/${i}.list ${StudyName}/${runPeriod}/${Channels}/${i}.root analysis_config.config
   echo ./ssb_analysis ${runPeriod}/${i}.list ${StudyName}/${runPeriod}/${Channels}/${i}.root config 

end
