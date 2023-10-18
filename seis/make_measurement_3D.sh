#!/bin/bash
#SBATCH --job-name=make_measurement_3
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=48
#SBATCH --time=02:00:00
#SBATCH --partition=normal

STARTTIME=$(date +%s)
echo "start time is : $(date +"%T")"

source ~/.bashrc
conda activate py3

export gen="python ../generate_path_files.py -p ../paths.yml -s ../settings.yml -e ../event_list"

events=`$gen list_events`
periods=`$gen list_period_bands`

export UCX_TLS="knem,dc_x"
nproc=480

for e in $events
do
    ibrun -n $nproc python make_measurement.py ${e} obsd_3 obsd_3D_crust
    echo ${e} Done........ >> events_done_3D
    ibrun -n $nproc python make_measurement.py ${e} obsd_25 obsd_glad
    echo ${e} Done........ >> events_done_glad
    ibrun -n $nproc python make_measurement.py ${e} obsd_1 obsd_1D_crust
    echo ${e} Done........ >> events_done_1D
done

ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")"
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."
