#!/bin/bash

# Validate input arguments
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <inputlist> <runPeriod> <StudyName> <Channels> <Sample> <SEDir>"
    exit 1
fi

# Assign input arguments
inputlists=($1)
runPeriod=$2
StudyName=$3
Channels=$4
Sample=$5
OutputSEDir=$6

echo $runPeriod
configdir=""
confch=""
echo "Good"

# Set config files w.r.t. Channels
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

# Create output directory if OutputSEDir is not "None"
if [ "$OutputSEDir" != "None" ]; then
    echo "Creating SE Directory: $OutputSEDir"
    mkdir -p $OutputSEDir/${StudyName}/${runPeriod}/${Channels}/${Sample}
    if [ $? -ne 0 ]; then
        echo "Failed to create SE Directory: $OutputSEDir"
        exit 1
    fi
else 
    echo "Creating Output Directory in AN Dir:"
    mkdir -p output/${StudyName}/${runPeriod}/${Channels}/${Sample}
     
fi

for i in "${inputlists[@]}"; do
   # command line
   echo ./ssb_analysis ${runPeriod}/${i}.list ${StudyName}/${runPeriod}/${Channels}/${i}.root ${configpath}${confch}
   ./ssb_analysis ${runPeriod}/${Sample}/${i}.list ${StudyName}/${runPeriod}/${Channels}/$Sample/${i}.root ${configpath}${confch} ${OutputSEDir}
done

