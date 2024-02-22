import os
import sys
import subprocess

cmd = 'qdel'
for i in range(2993839,2998634+1):
    stcmd = cmd + ' %s' %(i)
    os.system (stcmd)
print "finshed --"
