#!/bin/bash
#SBATCH --job-name=C6H4F2
#SBATCH --partition=LocalQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=C6H4F2_%j.out

export g16root=/opt/g16
export GAUSS_SCRDIR=$HOME/scratch
source $g16root/g16/bsd/g16.profile

g16 C6H4F2.com C6H4F2.log
