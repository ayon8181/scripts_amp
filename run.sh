#!/bin/bash
#SBATCH --job-name=workflow_test
#SBATCH --account=EAR21003
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --partition=small

source /home1/09038/ayon8181/.bashrc
conda activate py3

 python generate_path_files.py folders
 python generate_path_files.py converter
 cd converter
 ./convert_to_asdf.sh
 cd ..
 python generate_path_files.py proc
 cd proc
 ./run_preprocessing.sh
 cd ..
 python generate_path_files.py windows
 cd windows
 ./select_windows.sh
 cd ..
 python generate_path_files.py measure
 cd measure
 ./run_measureadj.sh
 cd ..
 python generate_path_files.py stations
 cd stations
 ./extract_stations.sh
 cd ..
 python generate_path_files.py filter
 cd filter
 ./filter_windows.sh
 cd ..
 python generate_path_files.py adjoint
 cd adjoint
 ./run_pyadj_mt.sh
 cd ..
 python generate_path_files.py weight_params
 python generate_path_files.py weight_paths
 cd weights
 ./calc_weights.sh
 cd ..
 python generate_path_files.py sum
 cd sum_adjoint
 ./sum_adjoint.sh
 cd ..
