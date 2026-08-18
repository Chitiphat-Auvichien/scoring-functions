#!/bin/bash
#SBATCH --job-name=C2H4
#SBATCH --partition=LocalQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=C2H4_%j.out

export g16root=/opt/g16
export GAUSS_SCRDIR=$HOME/scratch
source $g16root/g16/bsd/g16.profile

g16 C2H4.com C2H4.log
