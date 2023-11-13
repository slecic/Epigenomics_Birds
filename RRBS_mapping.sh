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

mkdir ${workdir}
mkdir ${outdir}

RUN=$1

#Trim adapters (for single-end reads)
trim_galore --fastqc -j 8 --quality 30 --rrbs --output_dir ${workdir} ${rawdata}/${RUN}*.fastq.gz

#Map reads
bismark --parallel 10 --output_dir ${workdir} --genome ${genomedir} ${workdir}/${RUN}_trimmed.fq.gz

#Extract methylation
bismark_methylation_extractor --parallel 6 --gzip --bedgraph --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${workdir}/ ${workdir}/${RUN}_trimmed_bismark_bt2.bam

#Only keep positions within predicted inserts
bedtools intersect -wb -a ${cutsites} -b ${workdir}/${RUN}*_trimmed_bismark_bt2.bismark.cov.gz | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort | uniq > ${workdir}/${RUN}.filt_cutsite.tmp

#Remove any positions overlapping a transition (C-T SNPs)
bedtools subtract -A -a ${workdir}/${RUN}.filt_cutsite.tmp -b ${c2tsnp} > ${workdir}/${RUN}.filt_cutsite_c2tsnp.tmp

#Only keep sites on the major chromosomes + and sort according to chromosome and position
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}.filt_cutsite_c2tsnp.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort -k1,1 -k2,2n | uniq > ${workdir}/${RUN}.filt_cutsite_c2tsnp_chrs.tmp

#Divide files further into repeat-region 5mC, and then our final 5mC positions for analysis (+ and sort according to chromosome and position)
bedtools intersect -wb -a ${repeats} -b ${workdir}/${RUN}.filt_cutsite_c2tsnp_chrs.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort -k1,1 -k2,2n | uniq | bgzip -c > ${outdir}/${RUN}.CpG_in_Repeats.cov.gz

#Subtract repeats for the FINAL file for our analysis
bedtools subtract -A -a ${workdir}/${RUN}.filt_cutsite_c2tsnp_chrs.tmp -b ${repeats} | bgzip -c > ${outdir}/${RUN}.CpG_5mC.cov.gz
