#!/bin/bash

source "/cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh"


# Validate input arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <dir> <inputlist> <runPeriod> <StudyName> <Channels> <Sample> <SEDir>"
    exit 1
fi

# Assign directory argument and change to that directory
dir=$1
shift # Remove the first argument so we can use the rest as before

# Determine the target directory
if [ "$dir" = "./" ] || [ "$dir" = "current" ]; then
    dir=$(pwd)
fi

# Change to the specified directory
cd "$dir" || { echo "Error: Cannot change to directory $dir"; exit 1; }

# Check if ssb_analysis exists in the specified directory
if [ ! -f "./ssb_analysis" ]; then
    echo "Error: ssb_analysis program not found in directory $dir"
    exit 1
fi

# Assign the rest of the input arguments
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
    mkdir -p "$OutputSEDir/${StudyName}/${runPeriod}/${Channels}/${Sample}"
    if [ $? -ne 0 ]; then
        echo "Failed to create SE Directory: $OutputSEDir"
        exit 1
    fi
else
    echo "Creating Output Directory in AN Dir:"
    mkdir -p "output/${StudyName}/${runPeriod}/${Channels}/${Sample}"
fi

for i in "${inputlists[@]}"; do
   # command line
   echo ./ssb_analysis "${runPeriod}/${i}.list" "${StudyName}/${runPeriod}/${Channels}/${i}.root" "${configpath}${confch}"
   ./ssb_analysis "${runPeriod}/${Sample}/${i}.list" "${StudyName}/${runPeriod}/${Channels}/$Sample/${i}.root" "${configpath}${confch}" "${OutputSEDir}"
done

