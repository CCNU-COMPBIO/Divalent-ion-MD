#!/bin/bash
#SBATCH --job-name=CA_20mM1
#SBATCH --account=rrg-panch
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -t 24:00:00                 # time (days-hours:minutes:seconds)

echo 'STARTING JOB' 	# prints to your output file

module load MistEnv/2021a cuda/10.2.2 gcc/8.5.0 anaconda3/2021.05 openmpi/4.0.5
export PATH=$HOME/amber22/bin:$PATH
export LD_LIBRARY_PATH=$HOME/amber22/lib:$LD_LIBRARY_PATH
export AMBERHOME=$HOME/amber22
source $AMBERHOME/amber.sh

NUM_PES=$(expr $SLURM_CPUS_PER_TASK - 1 )

source /home/p/panch/yunhuip/amber22/amber.sh

pmemd.cuda -O -i Min.in -o output/Min.out -p 1kx5_amber_1264_na_tip4pew.prmtop -c 1kx5_amber_1264_na_tip4pew.inpcrd -r output/Min.ncrst -inf output/Min.mdinfo
pmemd.cuda -O -i Equil_v.in -o output/Equil_v.out -p 1kx5_amber_1264_na_tip4pew.prmtop -c output/Min.ncrst -r output/Equil_v.ncrst -x output/Equil_v.nc -inf output/Equil_v.mdinfo
pmemd.cuda -O -i Equil_pt.in -o output/Equil_pt.out -p 1kx5_amber_1264_na_tip4pew.prmtop -c output/Equil_v.ncrst -r output/Equil_pt.ncrst -x output/Equil_pt.nc -inf output/Equil_pt.mdinfo
pmemd.cuda -O -i Prod.in -p 1kx5_amber_1264_na_tip4pew.prmtop -c output/Equil_pt.ncrst -o output/Prod.out -r output/Prod.ncrst -inf output/Prod.info -x output/Prod.nc

echo 'JOB ENDED'	# prints to your output file

#submit job
#sbatch --partition=gpu --gres=gpu:v100:1  --time=168:00:00  sub_job_GPU.sh
