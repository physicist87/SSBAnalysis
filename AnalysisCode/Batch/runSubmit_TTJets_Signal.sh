#!/bin/bash 
PBS -o /dev/null 
qsub -q cms  -l walltime=48:00:00,cput=48:00:00  -e /d3/scratch/sha/Analyses/SSB/MyAnalysis/25ns_Anlysis/Run2016/withMiniTree/Origin/Batch/Log/ -o /dev/null/Log/ -N TTJets_Signal_2 TTJets_Signal_2.sh
