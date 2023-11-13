#!/usr/bin/env Rscript


# Author:Sonja Lecic
# Originally written by Justin Meröndun
# Modified by Sonja 20.03.2021 to have the whole pipeline run faster 
# beacuse of a very large number of contings especially in case of the Blue tit genome

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
    # pattern that enzyme recognizes; example for HaeIII : GG'CC
    cs_5p1 = cs_5p1
    cs_3p1 = cs_3p1
    bluetit.dig <- insilico.digest(bluetit, cs_5p1, cs_3p1, verbose=T)
    bluetit.sel <- adapt.select(bluetit.dig, type="AA", cs_5p1, cs_3p1)
    nar.bluetit <- size.select(bluetit.sel,  min.size = minsize, max.size = maxsize, graph=F, verbose=T)
    dt <- as.data.frame(nar.bluetit@ranges)
    chr <- rep(gsub(pattern = "\\.SEQ.fasta$", "", (i)), nrow(dt))
    #chr <- rep(tools::file_path_sans_ext(i), nrow(dt))
    dtf <- cbind(dt, chr)
    bed <- dtf[, c(5, 2 ,3)]
    outputFile <- paste0(tools::file_path_sans_ext(i), ".bed")
    write.table(bed, file = outputFile, col.names = F, row.names = F, sep = "\t", quote = FALSE)
  }
}

## call the function with desired restriction pattern and fragment size
restrictionEnzyme(cs_5p1 = "C", cs_3p1 = "CGG", minsize = 200, maxsize = 400)


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
write.table(cutsites, file = "FragmentsMspI200_400.bed", col.names = F, row.names = F, sep = "\t", quote = FALSE)



