#!/bin/bash
#SBATCH --job-name=convert_to_asdf
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=03:00:00
#SBATCH --partition=small

STARTTIME=$(date +%s)
echo "start time is : $(date +"%T")"

source ~/.bashrc
conda activate py3
i=0
nproc=1
offset=0
offset_value=1
n_parallel_jobs=55
p=all

for f in paths/obsd_glad*.path.json
do
    ibrun -n $nproc -o $offset pypaw-convert_to_asdf -f $f -v -s &
    i=$((i+1))
    offset=$((offset+offset_value))
    if [[ $i -ge $n_parallel_jobs ]]; then
    wait
    i=0
    offset=0
   fi
done

ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")"
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."
