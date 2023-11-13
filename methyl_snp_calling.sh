#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=lecic@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

# The Files
# genome folder
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1
# Parus major genome fasta
genome=${genomedir}/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa
# putative cutsites simulated with simRad
cutsites=${genomedir}/CutsitesParusMspI100_400.bed
# list of all SNPs
snps=${genomedir}/Parus_SNPs.bed
# bed file with coordinates of the repeats
repeats=${genomedir}/Parus_Repeats.bed
# bed file with lengths of all major chromosomes
chrs=${genomedir}/Chromosomes.genome
# bed file with coordinates of major chromosomes
chrsbed=${genomedir}/Chromosomes.bed
# bed file with coordinates of the mitochondrion
mtdna=${genomedir}/mtDNA.bed
#bed file with 0-based coordinates for C-T and A-G SNPs
c2tsnp=${genomedir}/Parus_C2Tsnps.bed

#Directory paths
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0024/rawdata/Bisulfite-Seq/RRBS
workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/work
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/out

RUN=$1

# First we need to sort the bam files by reference position with picard
picard SortSam I=${workdir}/${RUN}_trimmed_bismark_bt2.bam O=${workdir}/${RUN}_trimmed_bismark_bt2_sort.bam SORT_ORDER=coordinate

# Then we will add read groups to sorted bam files
## 1) extract sample name and read group name
fbname=$(basename ${workdir}/${RUN}_trimmed_bismark_bt2_sort.bam | sed -E 's/(([^_]*_){4}).*/\1/; s/_$//' )
frgname=$(basename ${workdir}/${RUN}_trimmed_bismark_bt2_sort.bam | sed -E 's/(([^_]*_){6}).*/\2/; s/_$//' )
## 2) add read groups with picard
picard AddOrReplaceReadGroups I=${workdir}/${RUN}_trimmed_bismark_bt2_sort.bam O=${workdir}/${RUN}_trimmed_bismark_bt2_sort_RG.bam ID=`basename $fbname`  LB=`basename $frgname` PL=illumina PU=`basename $frgname` SM=`basename $fbname` CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate

# call SNPs with biscuit, compress and tabix index them
biscuit pileup -o ${outdir}/${RUN}_trimmed_bismark_bt2_sort_RG_biscuit.vcf ${genomedir}/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa ${workdir}/${RUN}_trimmed_bismark_bt2_sort_RG.bam

bgzip ${outdir}/${RUN}_trimmed_bismark_bt2_sort_RG_biscuit.vcf

tabix -p vcf ${outdir}/${RUN}_trimmed_bismark_bt2_sort_RG_biscuit.vcf.gz
