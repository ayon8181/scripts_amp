#!/bin/bash
#SBATCH --job-name=proc
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --time=02:00:00
#SBATCH --partition=small

STARTTIME=$(date +%s)
echo "start time is : $(date +"%T")"

source ~/.bashrc
conda activate py3

python plot_map.py


ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")"
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."