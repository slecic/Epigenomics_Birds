
##############################
### ANALYSIS IN methylKit #####
library(methylKit)

# create a list of all bismark cov files
file.list=list.files(path = "/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files", pattern = '[z]', full.names = F)
file.list=as.list(file.list)

# create a list of sample IDs from file names
sampleID=list()
for (i in 1:length(file.list)) {
  name = substring(file.list[[i]], 1, 16)
  print(name)
  sampleID=append(sampleID, name)
}

# read files in methylKit
myobj=methRead(file.list[c(1, 2, 3, 4, 13, 14, 15, 16, 17, 18, 19, 20, 5, 6, 7, 8, 10, 11, 12, 21, 22, 23, 24, 25, 26, 27, 28,
                           c(29:85))],
               pipeline = "bismarkCoverage",
               sample.id=sampleID[c(1, 2, 3, 4, 13, 14, 15, 16, 17, 18, 19, 20, 5, 6, 7, 8, 10, 11, 12, 21, 22, 23, 24, 25, 26, 27, 28, c(29:85))],
               assembly = "Parus.major.v1.1",
               treatment=c(rep(0, 12), rep(1, 72)))

length(myobj)

# filter for coverage
filtered.myobj=filterByCoverage(myobj,lo.count=20,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

# normalize coverage
filt.norm.myobj=normalizeCoverage(filtered.myobj, method="median")

# unite basic
#meth=methylKit::unite(filt.norm.myobj, destrand=FALSE, min.per.group=1L)
meth=methylKit::unite(filt.norm.myobj, destrand=FALSE)


### PCA #####
png("ParusandParusPublicandCya.png", width = 1500, height = 1000, units = "px", pointsize = 10)
PCASamples(meth)
dev.off()

methper <- percMethylation(meth)
methper[methper == 0] <- NA
methper <- na.omit(methper)
methpert <- t(methper)
pca.methpert <- prcomp(methpert, scale. = T)
pca.methpert_out <- as.data.frame(pca.methpert$x)
pca.methpert_out$Species <- c(rep("Blue tit", 12), rep("Great tit", 73))
pca.methpert_out$Environment <- c(rep("Germany", 23), rep("Netherlands", 62))
pca.methpert_out$Sex <- c(rep("male", 2), rep("female", 2), rep("male", 15), rep("female", 2), rep("male", 2), "female", "male", rep("unknown", 2), rep("female", 57))
pca.methpert_out$LifeStage <- c(rep("adult", 4), rep("chick", 2), rep("yearling", 2), rep("chick", 1), rep("yearling", 2), rep("chick", 2), rep("yearling", 2), rep("chick", 2), rep("yearling", 2), rep("adult", 6),
                                rep("chick", 2), rep("adult"), 36, rep("chick"), 4, rep("adult"), 17)
head(pca.methpert_out)

png("BluevsGreatTitplusPublicParus.png", width = 1200, height = 1000, units = "px", pointsize = 40)
autoplot(pca.methpert, data = pca.methpert_out, colour = 'Species', shape='LifeStage', 
         label =T, label.label = 'Sex', size = 8, label.vjust = 2, label.size = 6)+
  theme_light()+
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        legend.title=element_text(size=25),
        legend.text=element_text(size=25))+
  scale_color_manual(values = c("darkorange1","dodgerblue3"))
dev.off()



# calculate differentially methylated sites using Fisher's exact test and BH correction
myDiff=calculateDiffMeth(meth, adjust = "BH", test = "fast.fisher", slim = F)
# BH cutoff
myDiff001=getMethylDiff(myDiff,difference=25, qvalue = 0.05)
#cutoff <- -log10(myDiff001$qvalue)


### create a FUNCTION plotDMRs
data <- data.frame(chr = myDiff$chr,
                   win = myDiff$start,
                   pval= -log10(myDiff$qvalue),
                   mdiff = myDiff$meth.diff)

data=data[order(match(data$chr, chrlen$V1)),]
data[data == "Inf"] <- NA

chrlen <- read.table("/Volumes/LaCie/PhD/data/RRBSParus/cluster/download_files/Chromosomes.genome")
chrlen=chrlen[order(chrlen$V1),]
chromosomes=unique(data$chr)
chromosomes=as.character(chromosomes)
chromosomes_length=matrix(chrlen$V2, nrow = 1, dimnames = list(NULL, chromosomes))

plotDMPs <- function(data, chromosomes, chromosomes_length, chr_len, highlight=NULL, meth=NULL, demeth=NULL, manhatan = NULL){
  
  off <- 0
  upper_lim <- max(na.omit(data$mdiff))
  low_lim <- min(na.omit(data$mdiff))
  epimut <- c()
  epimutmeth <- c()
  epimutdemeth <- c()
  for(i in seq_along(chromosomes)){
    chr <- chromosomes[i]
    if(i%%2){col <- "deepskyblue4"}else{col <- "darkorange1"}
    temp <- data[which(data$chr == chr), ]
    print(head(temp))
    highlight = highlight
    snpsh <- highlight[which(highlight$chr == chr),]
    # table of DMPs or DMRs to highlihgt from the dataset
    mutationstb <- temp[which(temp$win %in% snpsh$win),]
    # save the number of snps for each chromosome for later
    mut <- nrow(mutationstb)
    
    # methylating epimutation
    meth = meth
    methchr <- meth[which(meth$chr == chr),]
    epimeth <- nrow(methchr)
    
    # demethylting epimuations
    demeth = demeth
    demethchr <- demeth[which(demeth$chr == chr),]
    epidemeth <- nrow(demethchr)
    
    if(chr == chromosomes[1]){
      plot(temp$win, temp$mdiff, xaxt = "n", xlim = c(0, sum(chromosomes_length[1, ])), ylim = c(upper_lim, low_lim), col = col, pch=19,
           ylab = "Change in % methylation", xlab = "Chromosomes", cex.lab=1.5)
    }else if(chr == chromosomes[1]){
      # highlight desired DMPs or DMRs
      points(mutationstb$win, mutationstb$pval,  col = "green4", pch = 19)
      # visualize regions with epi-methylation calls
      abline(v=methchr$midpos, col = "turquoise3", lty=1, lwd=1)
      # visualize regions with epi-demethylation calls
      abline(v=demethchr$midpos, col = "violetred3", lty=1, lwd=1)
      
    }else{
      points(temp$win+off, temp$mdiff, col = col, pch=19)
      points(mutationstb$win, mutationstb$pval,  col = "green4", pch = 19)
      abline(v=methchr$midpos+off, col = "turquoise3", lty=1, lwd=1)
      abline(v=demethchr$midpos+off, col = "violetred3", lty=1, lwd=1)
    }
    off <- off + chromosomes_length[1, chr]
    # fdr cutoff line
    abline(h=0, col = "grey")
    # label chromosomes on the x-axis
    chrl <- c(0, chrlen$V2[-32])
    v=c(0 + cumsum(chrl))
    axis(1, at=v[-length(v)] + diff(v) / 2, labels=chromosomes, las=2)
    manhatan = manhatan
    title(main = manhatan, cex.main=3)
    
    if (!is.null(highlight)) {
      # save the number of DMPs or DMRs for each chromosome
      epimut<-append(epimut, mut)
      # caclulate the total number of DMPs or DMRs genomewide
      epimutotal <- sum(epimut)
      print(epimutotal)
    }
    if (!is.null(meth)) {
      # save the number of epi-methylation calls for each chromosome
      epimutmeth<-append(epimutmeth, epimeth) 
      # caclulate the total number of epi-methylation calls genomewide
      epimutmethtot <- sum(epimutmeth)
      print(epimutmethtot)
    }
    if (!is.null(demeth)) {
      # save the number of epi-demethylation calls for each chromosome
      epimutdemeth <- append(epimutdemeth, epidemeth)
      # caclulate the total number of epi-demethylation calls genomewide
      epimuthdemethtot <- sum(epimutdemeth)
      print(epimuthdemethtot)
    }
  }
  # print a data frame with chromosomes and corresponding number of DMPs or DMRs
  if (!is.null(highlight)) {
    epimutations<-data.frame(chr = chromosomes, number = epimut)
    print(paste("DMPs/DMRs per chromosome:"))
    print(epimutations)
  }
  # print a data frame with chromosomes and corresponding number of epi-methylation calls
  if (!is.null(meth)) {
    epimethtations <- data.frame(chr = chromosomes, number = epimutmeth)
    print(paste("epi-methylation calls per chromosome:"))
    print(epimethtations)
  }
  # print a data frame with chromosomes and corresponding number of epi-demethylation calls
  if (!is.null(demeth)) {
    epidemethtations <- data.frame(chr = chromosomes, number = epimutdemeth)
    print(paste("epi-demethylation calls per chromosome:"))
    print(epidemethtations)
  }
}

png("ChickVsAdult2.png", width = 3000, height = 1200, units = "px", pointsize = 12)
plotDMPs(data = data, chromosomes = chromosomes[-32], chromosomes_length = chromosomes_length, chr_len = chrlen, 
         meth = NULL, demeth = NULL, highlight = NULL, manhatan = "Germany chick vs Germany yearling; base-pair resolution")
dev.off()


## get percent methylation
methper <- percMethylation(meth)
head(methper)
methper <- as.data.frame(methper)
head(methper)
dim(methper)

png("proba.png", width = 800, height = 500, units = "px", pointsize = 10)
plot(methper$C2R0242_BL_ADL_F[1:10], pch=19, col="salmon")
dev.off()

methper10 <- methper
mp10 <- data.frame(midparent = c(1:78975))
mp10$mother <- methper10$C2R0242_BL_ADL_F
mp10$father <- methper10$C2R0553_BL_ADL_M
mp10$midparent <- (methper10$C2R0242_BL_ADL_F + methper10$C2R0553_BL_ADL_M)/2
mp10$chick1 <- methper10$B2X6537_BL_CHK_M
mp10$yearling1 <- methper10$B2X6537_BL_YRL_M
mp10$chick2 <- methper10$B2X6540_BL_CHK_M
mp10$yearling2 <- methper10$B2X6540_BL_YRL_M


mpchick1 <- lm(chick1~midparent, data=mp10)
summary(mpchick1)

mpchick2 <- lm(chick2~midparent, data=mp10)


png("mpchick.png", width = 1000, height = 800, units = "px", pointsize = 6)
plot(mp10$midparent, mp10$chick1, pch=1, col="salmon")
points(mp10$midparent, mp10$chick2, pch=1, col="lightblue")
abline(mpchick1, col="darkred")
abline(mpchick2, col="darkblue")
dev.off()


momchick1 <- lm(chick1~mother, data=mp10)
momchick2 <- lm(chick2~mother, data=mp10)

png("motherchick.png", width = 1000, height = 800, units = "px", pointsize = 6)
plot(mp10$mother, mp10$chick1, pch=1, col="salmon")
points(mp10$mother, mp10$chick2, pch=1, col="lightblue")
abline(momchick1, col="darkred")
abline(momchick2, col="darkblue")
dev.off()

dadchick1 <- lm(chick1~father, data=mp10)
dadchick2 <- lm(chick2~father, data=mp10)

png("fatherchick.png", width = 1000, height = 800, units = "px", pointsize = 6)
plot(mp10$father, mp10$chick1, pch=1, col="salmon")
points(mp10$father, mp10$chick2, pch=1, col="lightblue")
abline(dadchick1, col="darkred")
abline(dadchick2, col="darkblue")
dev.off()

ep <- subset(samplefull, EPsire == 1)
table(ep$ID)
epwp <- subset(samplefull, EPsire == 0.5)
table(epwp$ID)

subset(samplefull, genfather=="B2X7195")

runinfo <- read.csv("/Volumes/LaCie/PhD/data/RRBSParusPublic/runinfo.csv")
head(runinfo, 10)
colnames(runinfo)
runsampname <- read.csv("/Volumes/LaCie/PhD/data/RRBSParusPublic/run_sample_name.csv")
head(runsampname)
