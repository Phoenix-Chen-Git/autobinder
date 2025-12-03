#!/bin/bash -l

#SBATCH -N 1
#SBATCH -p nitrogen
#SBATCH --gres=gpu:1
#SBATCH -o %x.log
#SBATCH -e %x.log

conda activate /xcfhome/yzmeng/miniconda3/envs/dl_binder_design

python /xcfhome/yzmeng/tools/design/dl_binder_design/af2_initial_guess/predict.py -pdbdir $1 -outpdbdir "$2"af2ig/ -scorefilename "$2"af2igpae
