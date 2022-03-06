#!/bin/bash
# parallel job using 1 node 20 cpu processors. 
#SBATCH -N 1
#SBATCH -n 20

#total time scheduled. If the job is running longer than this time limit, it will be terminated.
#SBATCH --time=7-00:00:00

#on CCAS, the longest time you can schedule is 7 days (7-00:00:00). See the CCAS webpage for more details:
#https://docs.mpcdf.mpg.de/doc/computing/clusters/systems/ExtraterrestrialPhysics/MPE-CCAS.html
#the type of notes that you are submitting to
#SBATCH -p ccas128

path=outputSMR3Dacc

srun athena -i athinput.binary3D -m 1 $path > $path/log
#srun athena -i athinput.binary3D -d $path > $path/log
