# Epigenomics
All things methylation

---
title: "Parus_major_RRBS_pipeline"
author: "Sonja Lecic"
date: "5/25/2021"
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

# Bisulfate Conversion

Parus major latest genome version (v1.1) from April 2020 can be found
[here](https://www.ncbi.nlm.nih.gov/assembly/GCF_001522545.3#/def_asm_Primary_Assembly).

Cyanistes carules latest genome version (v1) from January 2018 can be found [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_002901205.1/).

## File preparation
### Chromosomes

Extract chromosomes bed and genome files from *fa.fai. Chromosomes.bed file contains #chr #start #end. Chromosomes genome file contains #chr #FullLength.
```{}
# make a .fai file with samtools
samtools faidx Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa
# grep all the chromsomes and make a list
grep 'chr' Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.fai > Chromosomes.list
#Convert 1-based index to 0-based .bed format 
awk '{OFS="\t"}{print $1, "0", $2-1}' Chromosomes.list > Chromosomes.bed
#Create a genome length file for bedtools as well, called a 'genome' file 
awk '{OFS="\t"}{print $1, $2-1}' Chromosomes.list > Chromosomes.genome
```
### Mitochondrion
Create a .bed file from mitochondrion separately.
```{}
# grep the mitochondrion in a list
grep 'chrM' Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.fai > mtDNA.list
#Convert 1-based index to 0-based .bed format 
awk '{OFS="\t"}{print $1, "0", $2-1}' mtDNA.list > mtDNA.bed
```
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
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1

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
bioawk -c fastx '{print $name"\t0\t"length($seq)-1}' Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa > Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.bed
```
And, now the R script:
```{}
#!/usr/bin/env Rscript


# Author:Sonja Lecic
# Originally written by Justin Meröndun
# Modified by Sonja on 20.03.2021 to have the whole pipeline run faster 

#######################################
## In silico digest with SimRAD ###
#######################################
# set the working directory to a folder where all the files will be created and saved
setwd("/Volumes/LaCie/PhD/data/genomes/BlueTit/")

# call bash script from the working directory (find explanations of what the script does directly in the script)
x <- 'bash /Volumes/LaCie/PhD/data/RRBSParus/cluster/scripts/chrSelect.sh'
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
# give a name to your file and save it in your working directory
write.table(cutsites, file = "CutsitesParusMspI100_400.bed", col.names = F, row.names = F, sep = "\t", quote = FALSE)
```
..and the bash script *chrSelect.sh* that is called within the R script:
```{}
#!/bin/bash
### this script is selecting each scaffold from the genome fasta and saving it as a separate fasta file

for i in $(awk '{print $1}' /Volumes/LaCie/genome_assemblies_genome_fasta_Parus/ncbi-genomes-2021-05-20/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.bed); do
#simRAD can only hold on ~250 mb scaffolds at once, subselect scaffold
grep -w ${i} /Volumes/LaCie/genome_assemblies_genome_fasta_Parus/ncbi-genomes-2021-05-20/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.bed > ${i}.GRAB.bed
#use seqtk to grab this scaffold from the genome and save it as a separate file 
seqtk subseq /Volumes/LaCie/genome_assemblies_genome_fasta_Parus/ncbi-genomes-2021-05-20/Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa ${i}.GRAB.bed > ${i}.SEQ.fasta
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
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1

bismark_genome_preparation ${genomedir}
```

### Libraries list
Make a list of libraries for Parus major because I will make a script that will take each library in parallel from the raw directory and run it through the whole mapping pipeline.
```{}
ls *.fastq.gz | rev | cut -c10- | rev > /dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/Libraries.list
```

And, the genome folder looks like this now:
```{}
total 2.1G
drwxrws--- 4 ra57xaf pr53da-dss-0024 4.0K Jun  9 09:40 Bisulfite_Genome   # output of bismark genome preparation
-rw-rw---- 1 ra57xaf pr53da-dss-0024  557 May 26 13:35 Chromosomes.bed   # coordinates of major chromosomes
-rw-rw---- 1 ra57xaf pr53da-dss-0024  491 May 26 13:32 Chromosomes.genome   # total lengths of all major chromosomes
-rw-rw---- 1 ra57xaf pr53da-dss-0024 2.9M Jun  9 12:21 CutsitesParusMspI100_400.bed   # coodinates of the putative cutsites
-rw-rw---- 1 ra57xaf pr53da-dss-0024  816 Jun  9 11:43 Libraries.list   # list of Parus libraries
-rw-rw---- 1 ra57xaf pr53da-dss-0024   13 May 26 13:40 mtDNA.bed   # coodinates of mtDNA
-rw-rw---- 1 ra57xaf pr53da-dss-0024  14M Jun  9 13:08 Parus_C2Tsnps.bed   # coordinates of C-T SNPs
-rw-rw---- 1 ra57xaf pr53da-dss-0024 986M May 26 12:38 Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa # Parus genome fasta
-rw-rw---- 1 ra57xaf pr53da-dss-0024  34M Jun  9 11:35 Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.cat.gz
-rw-rw---- 1 ra57xaf pr53da-dss-0024 993M Jun  9 11:35 Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.masked
-rw-rw---- 1 ra57xaf pr53da-dss-0024 9.0M Jun  9 11:35 Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.out
-rw-rw---- 1 ra57xaf pr53da-dss-0024 5.9M Jun  9 11:35 Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.out.gff
-rw-rw---- 1 ra57xaf pr53da-dss-0024 2.5K Jun  9 11:35 Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa.tbl
-rw-rw---- 1 ra57xaf pr53da-dss-0024 2.7M Jun  9 11:52 Parus_Repeats.bed   # coodinates of the repeats
-rw-rw---- 1 ra57xaf pr53da-dss-0024 2.8M Jun  9 11:51 repeats_chrom.tmp
-rw-rw---- 1 ra57xaf pr53da-dss-0024 2.9M Jun  9 11:47 repeats.tmp1
drwxrws--- 2 ra57xaf pr53da-dss-0024 4.0K Jun  9 11:35 scripts_genome_prep   # folder with scripts
-rw-rw---- 1 ra57xaf pr53da-dss-0024    0 Jun  8 14:57 Snakefile
```

## Trimming, Mapping and Variant calling

### The full pipeline:
RRBS_Pipeline_CG.sh script
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
for sample in $(cat /dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/Libraries.list); do sbatch -J ${sample} RRBS_Pipeline_CG.sh ${sample}; done 
```

et voila! The final file .CpG_5mC.cov.gz looks like this:
```{}
The final file contains 6 columns: #Chromosome #Start #End #Methylation_Percentage #Count_Methylated #Count_Unmethylated

zcat C2R0242_BL_ADL_F__SRR10606832__SE__976ce27c8be6fb513c6153f49b03a3a1.CpG_5mC.cov.gz | head
chr01	4634	4634	88.4057971014493	183	24
chr01	4641	4641	88.4615384615385	184	24
chr01	4650	4650	82.6923076923077	172	36
chr01	4673	4673	97.1153846153846	202	6
chr01	4682	4682	72.7272727272727    8	 3
chr01	4711	4711	100	                1	 0
chr01	4803	4803	17.3913043478261	12	57
chr01	4817	4817	34.7826086956522	24	45
chr01	4860	4860	100	                1	 0
chr01	4865	4865	50	                1	 1
```

## CpG count summary

Concatenate all COUNT files into one to make an R input file:
```{}
awk 'FNR>1 || NR==1' *tCOUNTS* > COUNTS.input
```

And R code:
```{}
library(ggplot2)
library(ggsci)

counts <- read.table("/Volumes/LaCie/PhD/data/RRBSParus/cluster/COUNTS.input",header=T)
head(counts)

ggcounts <- data.frame(Filter = c(rep(colnames(counts[1]), nrow(counts)), 
                                  rep(colnames(counts[2]), nrow(counts)),
                                  rep(colnames(counts[3]), nrow(counts)),
                                  rep(colnames(counts[4]), nrow(counts)),
                                  rep(colnames(counts[5]), nrow(counts))),
                       Counts = c(unlist(counts[, c(1:5)], use.names = FALSE)),
                       Library = c(paste(rep(counts$LIBRARY_NAME))))

head(ggcounts)

png("ParusCpGcounts.png", width = 30, height = 30, units = "cm", res=600)
ggplot(ggcounts,aes(x= Library, fill= reorder(Filter,-Counts), y=Counts))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic(base_size=20)+
  labs(fill='Filter') +
  #scale_fill_manual(values=mypal)+
  scale_color_futurama()+
  ylab("Number of CpGs")+
  theme(axis.text.x=(element_text(angle=65,hjust=1)))
dev.off()

```


## Methylation bias along reads (M-bias)

Within /work/ there will be many files with an m-bias output showing methylation biases along read lengths, which could indicate systematic biases in the data.

The file contains outputs for CpGs, CHG and CHH for each sample and looks like this:
```{}
head B2X6537_BL_CHK_M__SRR10606831__SE__c0c278066ffc4182f5ab186575ff6d0d_trimmed_bismark_bt2.M-bias.txt
CpG context
===========
position	count methylated	count unmethylated	% methylation	coverage
1	5338097	8086304	39.76	13424401
2	634	812	43.85	1446
3	396	1030	27.77	1426
4	108026	335702	24.35	443728
5	123684	282038	30.48	405722
6	166299	425252	28.11	591551
7	176873	468122	27.42	644995
```

Combine all the outputs only for CpGs.
This command will merge all the files taking only CpGs, and remove the first 3 lines.
```{}
grep -A 53 'CpG' *M-bias* | sed '/CpG/d' | sed '/======/d' | sed '/--/d' | sed '/count/d' > ../out/M_BIAS_CpG.input
```
And R code:
```{}
mbias <- read.table("/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files/M_BIAS_CpG.input")

mbias$position <- paste0(substr(mbias$V1, 100,101))
mbias$V1 <- paste0(substr(mbias$V1, 1,29))
Mbiasnames <- c("Sample", "Count_M", "Count_U", "Percent_M", "Coverage", "Position")
names(mbias) <- Mbiasnames

# make a neat function to reorder columns (for later analysis also)
moveCol <- function(data, col2move, where = "last", beforeafter = NULL) {
  temp <- setdiff(names(data), col2move)
  x <- switch(
    where,
    first = data[c(col2move, temp)],
    last = data[c(temp, col2move)],
    before = {
      if (is.null(beforeafter)) stop("must specify beforeafter column")
      if (length(beforeafter) > 1) stop("beforeafter must be a single character string")
      data[append(temp, values = col2move, after = (match(beforeafter, temp)-1))]
    },
    after = {
      if (is.null(beforeafter)) stop("must specify beforeafter column")
      if (length(beforeafter) > 1) stop("beforeafter must be a single character string")
      data[append(temp, values = col2move, after = (match(beforeafter, temp)))]
    })
  x
}
mbiasdf = moveCol(data = mbias, col2move = "Position", where = "after", beforeafter = "Sample")

# change positions column to numeric
mbiasdf$Position <- as.numeric(mbiasdf$Position)

# plot M-bias
library(ggsci)
library(ggplot2)

ggplot(mbiasdf,aes(x=Position, col=Sample))+
  geom_line(aes(y=Percent_M),stat="identity",show.legend=FALSE, group = 12)+
  theme_classic(base_size=16)+
  scale_color_futurama()+
  coord_cartesian(ylim=c(0,100))

ggplot(mbiasdf,aes(x=Position,col=Sample))+
  geom_line(aes(y=Coverage),stat="identity",show.legend=FALSE, group = 12)+
  theme_classic(base_size=16)+
  scale_color_futurama()

```

## Analysis in methylKit

And R code:
```{}
library(methylKit)

# create a list of all bismark cov files
file.list=list.files(path = "/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files", pattern = '[gz]', full.names = F)
file.list=as.list(file.list)

# create a list of sample IDs from file names
sampleID=list()
for (i in 1:length(file.list)) {
  name = substring(file.list[[i]], 1, 29)
  print(name)
  sampleID=append(sampleID, name)
}

# read files in methylKit
myobj=methRead(file.list,
               pipeline = "bismarkCoverage",
               sample.id=sampleID,
               assembly = "Parus.major.v1.1",
               treatment=c(rep(1, 12)))
               

# filter for coverage above 20
filtered.myobj=filterByCoverage(myobj,lo.count=20,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)
                                
# normalize coverage
filt.norm.myobj=normalizeCoverage(filtered.myobj, method="median")

# unite all samples into one table
meth=methylKit::unite(filt.norm.myobj, destrand=FALSE, min.per.group = 9L)
head(meth)
dim(meth)
```

#### Differentially methylated positions (DMPs)
And R code:
```{}
# Calculate differentially methylated sites (not regions) using Fisher's exact test and BH correction
myDiff=calculateDiffMeth(meth, adjust = "BH", test = "fast.fisher", slim = F)
# BH cutoff
myDiff001=getMethylDiff(myDiff,difference=25,qvalue=0.01)
cutoff <- -log10(myDiff001$qvalue)


data <- data.frame(chr = myDiff$chr,
                   win = myDiff$start,
                   pval= -log10(myDiff$qvalue))

data=data[order(match(data$chr, chrlen$V1)),]
data[data == "Inf"] <- NA

chrlen <- read.table("/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files/Chromosomes.genome")
chrlen=chrlen[order(chrlen$V1),]
chromosomes=unique(data$chr)
chromosomes=as.character(chromosomes)
chromosomes_length=matrix(chrlen$V2, nrow = 1, dimnames = list(NULL, chromosomes))

snps=read.table("/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files/somatic.snps.list")
names(snps) <- c("chr", "win")

# FUNCTION to plot plotDMPs with an option to highlight specific snps
plotDMPs <- function(data, chromosomes, chromosomes_length, chr_len, highlight=NULL){
  
  off <- 0
  pval_lim <- max(na.omit(data$pval))
  epimut <- c()
  for(i in seq_along(chromosomes)){
    chr <- chromosomes[i]
    if(i%%2){col <- "deepskyblue4"}else{col <- "darkorange1"}
    temp <- data[which(data$chr == chr), ]
    print(head(temp))
    highlight = highlight
    snpsh <- highlight[which(highlight$chr == chr),]
    # mutations table
    mutationstb <- temp[which(temp$win %in% snpsh$win),]
    # save the number of mutations for each chromosome for later
    mut <- nrow(mutationstb)
    if(chr == chromosomes[1]){
      plot(temp$win, temp$pval, xaxt = "n", xlim = c(0, sum(chromosomes_length[1, ])), ylim = c(0, pval_lim), col = col, pch=19,
           ylab = "-log10(p)", xlab = "chromosomes", cex.lab=1.5)
    }else if(chr == chromosomes[1]){
      # highlight snps from the list
      points(mutationstb$win, mutationstb$pval,  col = "turquoise3", pch = 19)
      #points(temp$win[temp$win %in% snpsh$win], temp$pval[temp$win %in% snpsh$win], col = "turquoise3", pch = 19)
    }else{
      points(temp$win+off, temp$pval, col = col, pch=19)
      #points(temp$win[temp$win %in% snpsh$win]+off, temp$pval[temp$win %in% snpsh$win], col = "turquoise3", pch = 19)
      points(mutationstb$win+off, mutationstb$pval,  col = "turquoise3", pch = 19)
    }
    off <- off + chromosomes_length[1, chr]
    abline(h=min(cutoff), col = "red")
    chrl <- c(0, chrlen$V2)
    v=c(0 + cumsum(chrl))
    axis(1, at=v[-length(v)] + diff(v) / 2, labels=chromosomes, las=2)
    
    if (!is.null(highlight)) {
  # save the number of epimutations for each chromosome
    epimut<-append(epimut, mut)
  # caclulate the total number of epimutations genomewide
    epimutotal <- sum(epimut)
    print(epimutotal)
    }
  }
  if (!is.null(highlight)) {
  epimutations<-data.frame(chr = chromosomes, number = epimut)
  print(epimutations)
  }
}

png("B2X6537_minCpGperGroup1L.png", width = 3000, height = 1200, units = "px", pointsize = 12)
plotDMPs(data = data, chromosomes = chromosomes, chromosomes_length = chromosomes_length, chr_len = chrlen, highlight = posdemeth)
dev.off()


```

#### Diffrentially methylated regions (DMRs)
And R code:
```{}
# calculate differentially methylated sites using Fisher's exact test and BH correction
myDiff=calculateDiffMeth(meth, adjust = "BH", test = "fast.fisher", slim = F)
# BH cutoff
myDiff001=getMethylDiff(myDiff,difference=25,qvalue=0.01)
cutoff <- -log10(myDiff001$qvalue)

data <- data.frame(chr = myDiff$chr,
                   win = (myDiff$end - myDiff$start)/2 + myDiff$start,
                   pval= -log10(myDiff$pvalue))

data=data[order(match(data$chr, chrlen$V1)),]
data[data == "Inf"] <- NA

chrlen <- read.table("/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files/Chromosomes.genome")
chrlen=chrlen[order(chrlen$V1),]
chromosomes=unique(data$chr)
chromosomes=as.character(chromosomes)
chromosomes_length=matrix(chrlen$V2[-32], nrow = 1, dimnames = list(NULL, chromosomes))
off <- 0
pval_lim <- max(na.omit(data$pval))
png("proba.png", width = 3000, height = 1200, units = "px", pointsize = 12)
for(i in seq_along(chromosomes)){
  chr <- chromosomes[i]
  if(i%%2){col <- "deepskyblue4"}else{col <- "darkorange1"}
  temp <- data[which(data$chr == chr), ]
  print(head(temp))
  if(chr == chromosomes[1]){
    plot(temp$win, temp$pval, xaxt = "n", xlim = c(0, sum(chromosomes_length[1, ])), ylim = c(0, pval_lim), col = col, pch=19,
         ylab = "-log10(p)", xlab = "chromosomes", cex.lab=1.5)
  }else{
    points(temp$win+off, temp$pval, col = col, pch=20)
  }
  off <- off + chromosomes_length[1, chr]
  abline(h=min(cutoff), col = "red")
  chrl <- c(0, chrlen$V2[-32])
  v=c(0 + cumsum(chrl))
  axis(1, at=v[-length(v)] + diff(v) / 2, labels=chromosomes, las=2)
}
dev.off()

```

## Statistical analysis with lme4
```{}
library(lme4)
library(lmerTest)

### pepare methylation data frame for statistical analysis
head(meth)
dim(meth)

methdf <- methylKit::getData(meth)
samplelist <- meth@sample.ids

# make a function to get a table that combines methylation calls and trios information for statistical analysis
makeMethData <- function(methdf, samplelist, trios) {
  
  sampleName=vector()
  for (i in 1:length(samplelist)) {
    name = substring(samplelist[i], 1, 29)
    sampleName=append(sampleName, name)
  }
  
  logitTransform <- function(p) { log(p/(1-p)) }
  n <- 3
  cov <- methdf[5:ncol(methdf)][, c(TRUE, rep(FALSE, n - 1))]
  numC <- methdf[5:ncol(methdf)][, c(FALSE, TRUE, FALSE)]
  
  fC <- cbind(methdf[1:2],logitTransform(numC/cov))
  
  sampleName=as.character(sampleName)
  names(fC)[-c(1:2)] <- sampleName
  
  print(fC)
  
  # trios is a data set with trios ID and sex, petrnity, brood, nestbox, sampling year information
  # put everything together in a final data table
  epidata = data.frame(chr = c(rep(fC$chr, length(trios$offspring))),
                        pos = c(rep(fC$start, length(trios$offspring))),
                        fCoffspring = unlist(fC[, trios$offspring], use.names = FALSE),
                        fCmom = unlist(fC[, trios$mom], use.names = FALSE),
                        fCdad = unlist(fC[, trios$dad], use.names = FALSE),
                        sex = c(rep(trios$sex, each=nrow(fC))),
                        lifestage = c(rep(trios$lifestage, each=nrow(fC))),
                        paternity = c(rep(trios$paternity, each=nrow(fC))),
                        brood = c(rep(trios$brood, each=nrow(fC))),
                        nestbox = c(rep(trios$nestbox, each=nrow(fC))),
                        sampyear = c(rep(trios$sampyear, each=nrow(fC))))
  print(epidata)
  
}
# call the function
epistats = makeMethData(methdf = methdf, samplelist = samplelist, trios = trios)

# convert any Infinity or NaN values to NA 
epistats[epistats == "-Inf"] <- NA
epistats[epistats == "Inf"] <- NA
epistats[epistats == "NaN"] <- NA

# call an lmm model
epi.lmer <- lmer(fCoffspring ~  fCmom + (1|lifestage), data = epistats)
episum <- summary(epi.lmer)
# plot
plot(epi.lmer)
qqnorm(resid(epi.lmer))

# make a nice table out of summary results with stargazer library
library(stargazer)
# change the class to make it compatible for stargazer
class(epi.lmer) <- "lmerMod"
epi.lmertable <- stargazer(epi.lmer, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


```


# SNP calling
Call SNPs from the RRBS mapped files in order to get all SNPs (A, T, C and G).
Frist, I will sort the bam files and add read groups with picard. As there is no need to realign around indels (Bismark/Bowtie2 allows for gapped aligments). Then, base quality recalibration needs to be done as Bis-SNP heavily depends on base quality score for the genotype probability calculation, but Illumina sequencing reads’ raw base quality score could not accurately reflect the error rate of the base. I will then call SNPs with bis-snp.
RRBS_SNPcall.sh script:
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


```

And create a bash script for each sample:
```{}
for sample in $(cat /dss/dsslegfs01/pr53da/pr53da-dss-0024/projects/2021_EpigeneticsTits/Bisulfate_Conversion/Genomes/Parus.major.v1.1/Libraries.list); do sbatch -J ${sample} RRBS_SNPcall.sh ${sample}; done 
```

## Merging and Filtering
```{}

#first, merge all vcf files together
bcftools merge *.vcf.gz > parus.raw.snps.vcf.gz

# Check the total number of SNPs before filtering
zcat parus.raw.snps.vcf.gz | grep -v '#' | wc -l

# 10133064 SNPs
 
# Keep only sites with a maf 0.05, genotype min depth 10 and max depth 150, shared between 50% of samples
vcftools --gzvcf raw.snps.vcf.gz --recode --recode-INFO-all --bed ../Chromosomes.bed --maf 0.05 --minDP 10 --maxDP 150 --max-missing .5 --out dp10.maf05.mm5.snps

# Check the total number of SNPs after filtering
zcat dp10.maf05.mm5.filt.snps.vcf.gz| grep -v '#' |  wc -l


```

# Kinship matrix
Calculate kinship matrix with vcftools
```{}

vcftools --gzvcf dp10.maf05.mm5.filt.snps.vcf.gz --relatedness2

```

#### Somatic mutations
Here, wanted to try to call mutations between a chick sample and a yearling with bisquit pileup -S option
```{}
# here i sue -x opstion to specify mutations rate as 4*10-4 (for Arabidopsis thaliana)
biscuit pileup -x 0.0004 -S -o ../out/somatic_mode2.vcf ../Parus.major_genome_Abel_Illumina.HiSeq_AllPathsv.April2013_v1.1.fa  -T B2X6537_BL_YRL_M__SRR10606829__SE__c1629af0da0e72e8224c3dd0fe15817e_trimmed_bismark_bt2_sort_RG.bam -I B2X6537_BL_CHK_M__SRR10606849__SE__24c91df030d21a17e338f81b88c83423_trimmed_bismark_bt2_sort_RG.bam

# filtering
vcftools --gzvcf somatic_mode2.vcf.gz --recode --recode-INFO-all --bed ../Chromosomes.bed --maf 0.05 --minDP 10 --maxDP 150 --out somatic.mut.dp10.maf05.vcf

# get vcf with snp that have SS=2 in INFO column
cat somatic.mut.dp10.maf05.vcf.recode.vcf  | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' |grep 'SS=2' > somatic.snps

# save chromosome and position of those SNP in a text file
awk '{OFS="\t"}{print $1, $2}' somatic.snps > somatic.snps.list

```
I did not fine any mutations between chick and yearling samples! Mutations only happen in the germline??!!
