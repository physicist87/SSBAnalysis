import os
import sys
import subprocess

path = os.getcwd() 
path=path.replace("Batch", "")
os.getcwd()
workdir = os.getcwd()
print path
#print cddir
#os.chdir (path)
os.system ('pwd')
os.system ('mkdir -p Log')
os.system ('ls')
#runcmd = "./ssb_analysis DYJetsToLL_M_10To50_1.list DYJetsToLL_M_10To50_1.root 0"
#print runcmd
x = input("Num of Job ? : ")
SampleFile = raw_input(" Data or Sample Name ? : ")
SubFile = SampleFile + "_"
#qsubcmd = 'qsub -q short '
qsubcmd = 'qsub -q cms '

### making submit File ###
f1 = open( "runSubmit_%s.sh" %(SampleFile) , 'w')
f1.write('#!/bin/bash \n')
f1.write('PBS -o /dev/null \n')
nonstd = "/dev/null"
#f1.write('cd /d3/scratch/sha/Analyses/SSB/MyAnalysis/DYJetsToLL_M_10To50/dyjetstoll_m_50/ \n')
#cdpath = "cd " + path +
#f1.write('cd '+ path + " \n")

for i in range(1,x+1):
    runing = SampleFile + "_%s.sh" % (i)
    f = open( runing, 'w')
    f.write('#!/bin/bash \n')
    f.write('cd '+ path + " \n")
#    f.write('cd /d3/scratch/sha/Analyses/SSB/MyAnalysis/DYJetsToLL_M_10To50/dyjetstoll_m_50/ \n')
#    runcmd = "./ssb_analysis DYJetsToLL_M_10To50_%s.list DYJetsToLL_M_10To50_%s.root 0 " % (i,i)
    runcmd = "./ssb_analysis "+ SampleFile + "_%s.list " % (i) + SampleFile + "_%s.root " % (i) + " 0 "
    f.write(runcmd)
#    qsu = qsubcmd +"./Batch/" + runing + " -N " + SampleFile + "_%s " % (i)
#    qsu = qsubcmd + " -l walltime=72:00:00,cput=48:00:00 " + " -e " + workdir +"/Log/" + " -o " + workdir + "/Log/" + " -N " + SampleFile + "_%s " % (i) + runing 
    #qsu = qsubcmd + " -l walltime=72:00:00,cput=48:00:00 " + " -e " + workdir +"/Log/" + " -o " + nonstd + " -N " + SampleFile + "_%s " % (i) + runing 
#    qsu = qsubcmd + " -l walltime=00:45:00,cput=00:30:00 " + " -e " + workdir +"/Log/" + " -o " + nonstd + "/Log/" + " -N " + SampleFile + "_%s " % (i) + runing 
    qsu = qsubcmd + " -l walltime=48:00:00,cput=48:00:00 " + " -e " + workdir +"/Log/" + " -o " + nonstd + "/Log/" + " -N " + SampleFile + "_%s " % (i) + runing 
    print qsu
    f1.write(qsu+'\n')
    runchMod = "chmod 755 " + runing
    os.system (runchMod)
    f.close()
f1.close()
subchMod = "chmod 755 " + "runSubmit_%s.sh" % (SampleFile)
os.system (subchMod)
batchcmd = "./runSubmit_%s.sh" % (SampleFile)
print batchcmd
#os.system (batchcmd)
#os.system (runcmd)

