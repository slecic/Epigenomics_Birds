#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=lecic@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

#path to genome directory
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1

bismark_genome_preparation ${genomedir}
