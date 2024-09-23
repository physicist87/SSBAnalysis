#!/bin/bash

# Set up the environment
# Below is the ROOT Env setting at the KNU Tier3
source "/cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh"

# Define the main path and change to that directory
MainPath=$1
cd "${MainPath}"
echo "  [ Condor Exe Log ] Setted Main Work Directory"
echo "  [ Condor Exe Log ]   -->  ${MainPath}"

# Get the input arguments
inputListPath=$2
inputListName=$3
outputDir=$4
outputName=$5
outputPath=$6

# Create the output directory if it doesn't exist
# in ssb_analysis.cpp, you need to set the dsipcap path for using SE_Storage
mkdir -p "${outputPath}/${outputDir}"

# Run the analysis
echo "  [ Condor Exe Log ] Now jobs will be proceeded"
./ssb_analysis "./${inputListPath}/${inputListName}.list" "${outputDir}/${outputName}" 0
echo "  [ condor Exe Log ] Jobs done!!!"

