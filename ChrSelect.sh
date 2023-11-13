#!/bin/bash
### this script is selecting each scaffold from the genome fasta and saving it as a separate fasta file

for i in $(awk '{print $1}' /Volumes/LaCie/PhD/data/genomes/BlueTit/cyaCae2ChrOnly.bed); do
#simRAD can only hold on ~250 mb scaffolds at once, subselect scaffold
grep -w ${i} /Volumes/LaCie/PhD/data/genomes/BlueTit/cyaCae2ChrOnly.bed > ${i}.GRAB.bed
#use seqtk to grab this scaffold from the genome and save it as a separate file 
seqtk subseq /Volumes/LaCie/PhD/data/genomes/BlueTit/cyaCae2ChrOnly.fa ${i}.GRAB.bed > ${i}.SEQ.fasta
wait
# remove the GRAB.bed file from the folder
rm *.GRAB.bed
done
