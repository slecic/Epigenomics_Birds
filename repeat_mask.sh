#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=lecic@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

#specify path to the genome and genome directory
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1

RepeatMasker ${genome} -species chicken -nolow -pa 20 -gff -dir ${genomedir}
