# Epigenomics
All things methylation

---
title: "Cyanistes_RRBS"
author: "Sonja Lecic"
date: "6/28/2021"
output: 
   cleanrmd::html_document_clean:
     theme: minicss
     toc: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = F)
```


```{css}
.columns {display: flex;}
h1 {color: purple;}
h2 {color: darkorange;}
```

# Bisulfate conversion

## Genome annotation

Cyanistes caerules latest genome version (v1) from January 2018 can be found [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_002901205.1/).

## File preparation
### Chromosomes

Extract chromosomes bed and genome files from *fa.fai. Chromosomes.bed file contains #chr #start #end. Chromosomes genome file contains #chr #FullLength.
```{}
# make a .fai file with samtools
samtools faidx Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa
# grep all the chromsomes and make a list
grep 'chr' Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.fai > Chromosomes.list
#Convert 1-based index to 0-based .bed format 
awk '{OFS="\t"}{print $1, "0", $2-1}' Chromosomes.list > Chromosomes.bed
#Create a genome length file for bedtools as well, called a 'genome' file 
awk '{OFS="\t"}{print $1, $2-1}' Chromosomes.list > Chromosomes.genome
```


### Mitochondria
There is no Mitochondrion characterized for Cynistes. I will map to Parus mtDNA when the time comes.

### Repeats
Create a bed file containing repeats. These will be analyzed separately. To extract the repeat I will use RepeatMasker.
```{}
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=lecic@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

#specify path to the genome and genome directory
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1/Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1

RepeatMasker ${genome} -species chicken -nolow -pa 20 -gff -dir ${genomedir}
```
Convert it to a standard bed:
```{}
awk '{print $1, $4, $5, $7, $10}' Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.out.gff > repeats.tmp1
grep -E -- 'chr' repeats.tmp1 > repeats_chrom.tmp
awk '{OFS="\t"}{print $1, $2, $3, $4, $5}' repeats_chrom.tmp | sed 's/"//g' > Parus_Repeats.bed
```

### Cut-sites
Identify cut-sites in the genome. For this, I will use a combination of an R and bash scripts.The script uses simRAD R package to find the cut-sites. I we are going to use the *MspI restriction enzyme* and *fragment size selection of 100-400bp* this needs to be put in the function.
Cut-sites bed file will be used to after mapping of the RRBS reads to subset the reads that overlap with these fragments in order to reduce the number erroneously mapped reads.

But, first we need to prepare the bed file for the whole genome (chromosomes + scaffolds; we can ectract from the final file with cutsites whatever we want then). For this I'll use bioawk to make a bed file (#chrom #start #end) in a 0 based bed format:
```{}
bioawk -c fastx '{print $name"\t0\t"length($seq)-1}' Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa > Genome.bed
```
And, now the RestEnzymesRRBS.R script:
```{}
#!/usr/bin/env Rscript


# Author:Sonja Lecic
# Originally written by Justin Meröndun
# Modified by Sonja 20.03.2021 to have the whole pipeline run faster 

#######################################
## In silico digest with SimRAD ###
#######################################
# set the working directory to a folder where all the files will be created and saved
setwd("/Volumes/LaCie/PhD/data/RRBSCyanistes/genome")

# call bash script from the working directory (find explanations of what the script does directly in the script)
x <- 'bash /Volumes/LaCie/PhD/data/RRBSCyanistes/chrSelect.sh'
system(x)

# load SimRAD library
library(SimRAD)
# make a function to create text a file with fragment ranges that takes as input: 
# - 5' and 3' nucleotide pattern that restriction enzyme recognizes
# - min and max size of selection of the fragments
restrictionEnzyme <- function(cs_5p1, cs_3p1, minsize, maxsize) {
  # make a list of all SEQ.fasta files
  myFiles <- list.files(pattern = "\\.SEQ.fasta$") 
  for(i in myFiles){
    bluetit <- ref.DNAseq(i, subselect.contigs=FALSE)
    # pattern that enzyme recognizes; for MspI : C'CGG
    cs_5p1 = cs_5p1
    cs_3p1 = cs_3p1
    bluetit.dig <- insilico.digest(bluetit, cs_5p1, cs_3p1, verbose=T)
    bluetit.sel <- adapt.select(bluetit.dig, type="AA", cs_5p1, cs_3p1)
    nar.bluetit <- size.select(bluetit.sel,  min.size = minsize, max.size = maxsize, graph=F, verbose=T)
    dt <- as.data.frame(nar.bluetit@ranges)
    chr <- rep(gsub(pattern = "\\.SEQ.fasta$", "", (i)), nrow(dt))
    dtf <- cbind(dt, chr)
    bed <- dtf[, c(5, 2 ,3)]
    outputFile <- paste0(tools::file_path_sans_ext(i), ".bed")
    write.table(bed, file = outputFile, col.names = F, row.names = F, sep = "\t", quote = FALSE)
  }
}

## call the function with desired restriction pattern and fragment size
restrictionEnzyme(cs_5p1 = "C", cs_3p1 = "CGG", minsize = 100, maxsize = 400)


## remove empty files if there are any (because some parts of a chromosome will have no cutsites and therefore produce emty bed files)
docs <- list.files(pattern = "\\.SEQ.bed$")   
inds <- file.size(docs) == 0 
file.remove(docs[inds])

# concatenate all the beds into one bed file with #Chrom #Start #End as columns
myBeds <- list.files(pattern = "\\.SEQ.bed$")
cutsites <-  read.table(myBeds[1], header = F)
for (f in myBeds[-1]){
  cutsites <- rbind(cutsites, read.table(f))
}
# give a name to your cutsites file and save it in your working directory
write.table(cutsites, file = "CutsitesCyanistesMspI100_400.bed", col.names = F, row.names = F, sep = "\t", quote = FALSE)



```
..and the bash script *chrSelect.sh* that is called within the R script:
```{}
#!/bin/bash                                                                                                                                 
### this script is selecting each scaffold from the genome fasta and saving it as a separate fasta file

for i in $(awk '{print $1}' /Volumes/LaCie/PhD/data/RRBSCyanistes/genome/Chromosome.bed); do
#simRAD can only hold on ~250 mb scaffolds at once, subselect scaffold                                                                      
grep -w ${i} /Volumes/LaCie/PhD/data/RRBSCyanistes/genome/Chromosome.bed > ${i}.GRAB.bed
#use seqtk to grab this scaffold from the genome and save it as a separate file                                                             
seqtk subseq /Volumes/LaCie/PhD/data/RRBSCyanistes/genome/Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY201\
4_v1.fa ${i}.GRAB.bed > ${i}.SEQ.fasta
wait
# remove the GRAB.bed file from the folder                                                                                                  
rm *.GRAB.bed

done

```

## Genome preparation

Genome will be mapped with bismark. Before the actuall mapping the genome needs to be prepared:
```{}
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
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1

bismark_genome_preparation ${genomedir}
```

### Libraries list
Make a list of libraries for Parus major because I will make a script that will take each library in parallel from the raw directory and run it through the whole mapping pipeline.
```{}
ls *.fastq.gz | rev | cut -c10- | rev > /dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1/Libraries.list
```
And, the genome folder looks like this now:
```{}
total 2.4G
drwxrws--- 4 ra57xaf pr53da-dss-0024 4.0K Jun 28 13:02 Bisulfite_Genome
-rw-rw---- 1 ra57xaf pr53da-dss-0024 8.7K Jun 28 12:07 Chromosomes.bed
-rw-rw---- 1 ra57xaf pr53da-dss-0024 8.0K Jun 28 12:07 Chromosomes.genome
-rw-rw---- 1 ra57xaf pr53da-dss-0024  15K Jun 28 12:07 Chromosomes.list
-rw-rw---- 1 ra57xaf pr53da-dss-0024 3.0M Jun 29 11:22 CutsitesCyanistesMspI100_400.bed # cutsites file
-rw-rw---- 1 ra57xaf pr53da-dss-0024 1.2G Jun 28 11:59 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa
-rw-rw---- 1 ra57xaf pr53da-dss-0024  37M Jun 28 12:45 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.cat.gz
-rw-rw---- 1 ra57xaf pr53da-dss-0024 1.1M Jun 28 12:06 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.fai
-rw-rw---- 1 ra57xaf pr53da-dss-0024 1.2G Jun 28 12:45 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.masked
-rw-rw---- 1 ra57xaf pr53da-dss-0024  10M Jun 28 12:45 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.out
-rw-rw---- 1 ra57xaf pr53da-dss-0024 6.9M Jun 28 12:45 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.out.gff
-rw-rw---- 1 ra57xaf pr53da-dss-0024 2.5K Jun 28 12:45 Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa.tbl
-rw-rw---- 1 ra57xaf pr53da-dss-0024 3.3M Jun 28 12:57 Cyanistes_Repeats.bed
-rw-rw---- 1 ra57xaf pr53da-dss-0024 632K Jun 28 12:12 Genome.bed 
-rw-rw---- 1 ra57xaf pr53da-dss-0024 1.9K Jun 28 13:05 Libraries.list # list of sample names
-rw-rw---- 1 ra57xaf pr53da-dss-0024  14M Jun 28 12:09 Parus_C2Tsnps.bed  # C to T SNP list
-rw-rw---- 1 ra57xaf pr53da-dss-0024 3.4M Jun 28 12:57 repeats_chrom.tmp
-rw-rw---- 1 ra57xaf pr53da-dss-0024 3.6M Jun 28 12:57 repeats.tmp1
drwxrws--- 2 ra57xaf pr53da-dss-0024 4.0K Jun 28 13:02 scripts
```

## Trimming, Mapping and Variant calling

### The full pipeline:
RRBS_CGpipeline_Cyanistes.sh script
```{}
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
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1
# Parus major genome fasta
genome=${genomedir}/Cyanistes.caeruleus_genome_BT333.1_Illumina.HiSeq_CeleraAssemblerv7.IDBAUDvMAY2014_v1.fa
# putative cutsites simulated with simRad
cutsites=${genomedir}/CutsitesCyanistesMspI100_400.bed
# bed file with coordinates of the repeats
repeats=${genomedir}/Cyanistes_Repeats.bed
# bed file with lengths of all major chromosomes
chrs=${genomedir}/Chromosomes.genome
# bed file with coordinates of major chromosomes
chrsbed=${genomedir}/Chromosomes.bed
# bed file with coordinates of the mitochondrion
mtdna=${genomedir}/mtDNA_Parus.bed
#bed file with 0-based coordinates for C-T and A-G SNPs
c2tsnp=${genomedir}/Parus_C2Tsnps.bed

#Directory paths
rawdata=/dss/dsslegfs01/pr53da/pr53da-dss-0024/rawdata/Bisulfite-Seq/RRBS
workdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1/work
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1/out

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

###########

### SUMMARIZE COUNTS
#Starting, raw positions
echo "RAW" >> ${outdir}/${RUN}.COUNTS.txt
zcat ${workdir}/${RUN}*_trimmed_bismark_bt2.bismark.cov.gz | wc -l >> ${outdir}/${RUN}.COUNTS.txt

#Positions remaining after cut-site filter
echo "CUTSITE_FILTER" >> ${outdir}/${RUN}.COUNTS.txt
cat ${workdir}/${RUN}.filt_cutsite.tmp | wc -l >> ${outdir}/${RUN}.COUNTS.txt

#Positions remaining after C-T SNP filter
echo "CT_SNP_FILTER" >> ${outdir}/${RUN}.COUNTS.txt
cat ${workdir}/${RUN}.filt_cutsite_c2tsnp.tmp | wc -l >> ${outdir}/${RUN}.COUNTS.txt

#Positions remaining on major chromosomes
echo "CHROMOSOME_FILTER" >> ${outdir}/${RUN}.COUNTS.txt
cat ${workdir}/${RUN}.filt_cutsite_c2tsnp_chrs.tmp | wc -l >> ${outdir}/${RUN}.COUNTS.txt

#Positions remaining after filtering sites within repeats
echo "REPEAT_FILTER" >> ${outdir}/${RUN}.COUNTS.txt
zcat ${outdir}/${RUN}.CpG_5mC.cov.gz | wc -l >> ${outdir}/${RUN}.COUNTS.txt

#Add a name column for each library for easier manipulation later
echo "LIBRARY_NAME" >> ${outdir}/${RUN}.COUNTS.txt
fbname=$(basename ${outdir}/${RUN}.COUNTS.txt | sed -E 's/(([^_]*_){6}).*/\1/; s/_$//' )
echo "$fbname" >> ${outdir}/${RUN}.COUNTS.txt

#Transpose the table so that we have 6 columns and one row with counts
cat ${outdir}/${RUN}.COUNTS.txt | pr -ts' ' --column 6 > ${outdir}/${RUN}.tCOUNTS.txt

```

And create a bash script for each sample:
```{}
#!/bin/bash

for sample in $(cat /dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Cyanistes.caeruleus.v1/Libraries.list); do sbatch -J ${sample} RRBS_CGpipeline_Cyanistes.sh ${sample}; done 
```
