#!/bin/bash
#SBATCH --job-name=make_measurement
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=48
#SBATCH --time=02:0:00
#SBATCH --partition=development

STARTTIME=$(date +%s)
echo "start time is : $(date +"%T")"

source ~/.bashrc
conda activate py3

export gen="python ../generate_path_files.py -p ../paths.yml -s ../settings.yml -e ../event_list"

events=`$gen list_events`
periods=`$gen list_period_bands`

export UCX_TLS="knem,dc_x"
nproc=384

for e in $events
do
    echo $e >> 1D_crust
    #ibrun -n $nproc python make_measurement_2.py ${e} obsd_3 obsd_3D_crust 0
    #ibrun -n $nproc python make_measurement_2.py ${e} obsd_25 obsd_glad 0
    #ibrun -n $nproc python make_measurement_2.py ${e} prem_3 prem_3D_crust 0
    #ibrun -n $nproc python make_measurement.py ${e} obsd_3 obsd_3D_crust 0
    #ibrun -n $nproc python make_measurement.py ${e} obsd_25 obsd_glad 0
    #ibrun -n $nproc python make_measurement.py ${e} prem_3 prem_3D_crust 0
    #ibrun -n $nproc python make_measurement_2.py ${e} prem_33 prem_3D_atten 0
    #ibrun -n $nproc python make_measurement.py ${e} ref_1 1D_ref 0
    #ibrun -n $nproc python make_measurement_2.py ${e} obsd_1 obsd_1D_crust 0
    #ibrun -n $nproc python make_measurement_2.py ${e} obsd_32 obsd_3D_crust_2 0
    #ibrun -n $nproc python make_measurement_2.py ${e} obsd_315 obsd_3D_crust_15 0
    #ibrun -n $nproc python make_measurement_2.py ${e} real_data real_data 0
    ibrun -n $nproc python make_measurement_2.py ${e} obsd_31 obsd_3D_crust_1 0
    ibrun -n $nproc python make_measurement_2.py ${e} obsd_318 obsd_3D_crust_18 0
    ibrun -n $nproc python make_measurement.py ${e} obsd_31 obsd_3D_crust_1 0
    ibrun -n $nproc python make_measurement.py ${e} obsd_318 obsd_3D_crust_18 0
    #ibrun -n $nproc python make_measurement_2.py ${e} prem_16 prem_q16 0
done

ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")"
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."
