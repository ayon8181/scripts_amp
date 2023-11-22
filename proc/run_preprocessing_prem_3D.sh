#!/bin/bash
#SBATCH --job-name=proc
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=56
#SBATCH --time=00:20:00
#SBATCH --partition=normal

STARTTIME=$(date +%s)
echo "start time is : $(date +"%T")"

source ~/.bashrc
conda activate py3

export gen="python ../generate_path_files.py -p ../paths.yml -s ../settings.yml -e ../event_list"

events=`$gen list_events`
periods=`$gen list_period_bands`
nproc=56
nworkers=$SLURM_NNODES
jobfile=jobs_proc_prem_3D_crust_${SLURM_JOBID}
rm -rf ${jobfile}
touch $jobfile


export UCX_TLS="knem,dc_x"


for e in $events
do
    for p in $periods
       do
#	if [[ -f `$gen filename proc_obsd_glad $e $p` ]]; then
#	    echo "skip:" $e $p
#	    continue
#	fi
	echo running: $i $e $p
 	echo python proc.py \
    	   -p ./parfile/proc_prem_3D_crust.${p}.param.yml \
    	   -f ./paths/proc_prem_3D_crust.${e}.${p}.path.json >> ${jobfile}
    	  
	


    done
done
../run_mpi_queue.py $nproc $nworkers ${jobfile}

ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")"
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."
