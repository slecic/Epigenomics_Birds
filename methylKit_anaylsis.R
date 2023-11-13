
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
