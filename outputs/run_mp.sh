#!/bin/bash
#SBATCH --job-name=proc
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --partition=development
#SBATCH --output=proc.out

STARTTIME=$(date +%s)
echo "start time is : $(date +"%T")"

source ~/.bashrc
rm -rf /home1/09038/ayon8181/.gmt/sessions
conda activate pygmt
export UCX_TLS="knem,dc_x"
#python event_wise_map.py
#python plot_real_mp.py
#python plot_relative.py
#python plot_seis_2.py
python convert_to_phase_real.py
#python convert_to_phase.py
#python plot_avg.py

ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")"
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."
