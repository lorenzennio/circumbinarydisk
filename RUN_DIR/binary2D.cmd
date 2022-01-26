#!/bin/bash
# parallel job using 1 node 20 cpu processors. 
#SBATCH -N 6
#SBATCH -n 112
#total time scheduled. If the job is running longer than this time limit, it will be terminated.
#on CCAS, the longest time you can schedule is 7 days (7-00:00:00). See the CCAS webpage for more details:
#https://docs.mpcdf.mpg.de/doc/computing/clusters/systems/ExtraterrestrialPhysics/MPE-CCAS.html
#the type of notes that you are submitting to
#SBATCH -p ccas128

srun athena -i athinput.binary2D -m 1 outputSMR2D/ > outputSMR2D/log
srun athena -i athinput.binary2D -d outputSMR2D/ > outputSMR2D/log
