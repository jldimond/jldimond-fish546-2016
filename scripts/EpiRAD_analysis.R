##################################################################
## This script reads in a text file derived from the "Base Counts"
## .vcf output from ipyrad. Base counts are used for analysis of 
## EpiRADseq data.

# Read in data file
data <- read.delim("data3-2.txt", header=FALSE)
# Since the base counts were split into four columns, these need 
# to be summed
data2 <- t(sapply(seq(4,ncol(data), by=4), function(i) {
     indx <- i:(i+3)
     rowSums(data[indx[indx <= ncol(data)]])}))
#The resulting file needs to be transposed and turned into a dataframe
data3 <- as.data.frame(t(data2))
#Add header names
header <- read.delim("header_data3.txt", header=FALSE)
names <- as.vector(t(header))
colnames(data3) <- names
#Select samples of interest (some have very low sample sizes)
data4 <- data3[,c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
                  22,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,
                  41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56)]

#Remove ddr rows that have any zeros. The premise here is that zeros 
#in the EpiRAD dataset are informative because they may reflect 
#methylation, but they could also relfect true absence of the locus
#in the library. Here the ddRAD library serves to standarize the EpiRAD
#library. Any zeros in the ddRAD libary are treated as absence of the
#locus, thereby leaving zeros in the EpiRAD library only where the 
#locus was counted in the ddRAD library.

data5 <- data4[apply(data4[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,
                             31,33,35,37,39,41,43,45,47,49,51)],1,
                     function(z) !any(z==0)),] 


#We may want to look at a matrix with only EpiRAd data
data6 <- data5[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,
                  38,40,42,44,46,48,50,52)]

#################################################################
# Now use edgeR package to standardize count data by library size

library("edgeR")
#read in the file to edgeR
counts <- DGEList(counts=data5)
counts$samples
#TMM normalization (corrects for library size)
counts2 <- calcNormFactors(counts)
counts2$samples
#extract normalized counts
counts2_cpm <- cpm(counts2, normalized.lib.sizes=TRUE, log=TRUE)

##Plots to show ddRAD vs EpiRAD library (before normalization)
par(mfrow = c(7, 4))
par(mar = c(2,2,2,2)) 
par(oma = c(2,2,1,1)) 

for (i in seq(1,51, by = 2)){
  plot(data5[,i], data5[,i+1], xlab= colnames(data5[,i]), 
  ylab= colnames(data5[,i+1]))
}


#plot normalized counts
par(mfrow = c(7, 4))
par(mar = c(2,2,2,2)) 
par(oma = c(2,2,1,1)) 

for (i in seq(1,51, by = 2)){
  plot(counts2_cpm[,i], counts2_cpm[,i+1], xlab= colnames(counts2_cpm[,i]), 
       ylab= colnames(counts2_cpm[,i+1]))
}

##################################################################
#Using lm to get residuals

lm103 <- lm(counts2_cpm[,2] ~ counts2_cpm[,1])
lm104 <- lm(counts2_cpm[,4] ~ counts2_cpm[,3])
lm105 <- lm(counts2_cpm[,6] ~ counts2_cpm[,5])
lm106 <- lm(counts2_cpm[,8] ~ counts2_cpm[,7])
lm107 <- lm(counts2_cpm[,10] ~ counts2_cpm[,9])
lm108 <- lm(counts2_cpm[,12] ~ counts2_cpm[,11])
lm109 <- lm(counts2_cpm[,14] ~ counts2_cpm[,13])
lm110 <- lm(counts2_cpm[,16] ~ counts2_cpm[,15])
lm111 <- lm(counts2_cpm[,18] ~ counts2_cpm[,17])
lm112 <- lm(counts2_cpm[,20] ~ counts2_cpm[,19])
lm114 <- lm(counts2_cpm[,22] ~ counts2_cpm[,21])
lm115 <- lm(counts2_cpm[,24] ~ counts2_cpm[,23])
lm116 <- lm(counts2_cpm[,26] ~ counts2_cpm[,25])
lm117 <- lm(counts2_cpm[,28] ~ counts2_cpm[,27])
lm118 <- lm(counts2_cpm[,30] ~ counts2_cpm[,29])
lm121 <- lm(counts2_cpm[,32] ~ counts2_cpm[,31])
lm122 <- lm(counts2_cpm[,34] ~ counts2_cpm[,33])
lm123 <- lm(counts2_cpm[,36] ~ counts2_cpm[,35])
lm124 <- lm(counts2_cpm[,38] ~ counts2_cpm[,37])
lm125 <- lm(counts2_cpm[,40] ~ counts2_cpm[,39])
lm126 <- lm(counts2_cpm[,42] ~ counts2_cpm[,41])
lm127 <- lm(counts2_cpm[,44] ~ counts2_cpm[,43])
lm128 <- lm(counts2_cpm[,46] ~ counts2_cpm[,45])
lm129 <- lm(counts2_cpm[,48] ~ counts2_cpm[,47])
lm130 <- lm(counts2_cpm[,50] ~ counts2_cpm[,49])
lm131 <- lm(counts2_cpm[,52] ~ counts2_cpm[,51])

resid103 <- residuals(lm103)
resid104 <- residuals(lm104)
resid105 <- residuals(lm105)
resid106 <- residuals(lm106)
resid107 <- residuals(lm107)
resid108 <- residuals(lm108)
resid109 <- residuals(lm109)
resid110 <- residuals(lm110)
resid111 <- residuals(lm111)
resid112 <- residuals(lm112)
resid114 <- residuals(lm114)
resid115 <- residuals(lm115)
resid116 <- residuals(lm116)
resid117 <- residuals(lm117)
resid118 <- residuals(lm118)
resid121 <- residuals(lm121)
resid122 <- residuals(lm122)
resid123 <- residuals(lm123)
resid124 <- residuals(lm124)
resid125 <- residuals(lm125)
resid126 <- residuals(lm126)
resid127 <- residuals(lm127)
resid128 <- residuals(lm128)
resid129 <- residuals(lm129)
resid130 <- residuals(lm130)
resid131 <- residuals(lm131)


resid_all <- cbind(resid103, resid104, resid105, resid106, resid107, resid108,
                   resid109, resid110, resid111, resid112, resid114, resid115,
                   resid116, resid117, resid118, resid121, resid122, resid123,
                   resid124, resid125, resid126, resid127, resid128, resid129,
                   resid130, resid131)

#plot residuals
par(mfrow = c(7, 4))
par(mar = c(2,2,2,2)) 
par(oma = c(2,2,1,1)) 

for (i in 1:26){
  plot(resid_all[,i])
}

########################################################
#Read in sample info
sinfo <- read.table("sample_info.txt", colClasses = 'character', header = FALSE)
#without sample 101
sinfo2 <- sinfo[2:27,]
#transpose
tsinfo <- t(sinfo2)
#create vectors for morphotype and diameter (note whether sample 101
#was included or not)
type <- tsinfo[3,]
diam <- tsinfo[5,]
#without sample 112
type2 <- tsinfo[3,c(1:9,11:26)]
diam2 <- tsinfo[5,c(1:9,11:26)]
depth <- tsinfo[2,c(1:9,11:26)]

########################################################
#MDS
resid_t <- t(resid_all)
d <- dist(resid_t) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y)
text(x, y, labels = row.names(resid_t), col = rainbow(length(diam))[rank(diam2)], cex=.7)

nmds <- isoMDS(d)
plot(nmds$points, type = "n")
text(nmds$points, labels = row.names(resid_t))


########################################################
#PCA
pca <- prcomp(resid_t)
summary(pca)
eig <- (pca$sdev)^2

#plot
biplot(pca)

# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- pca$rotation
sdev <- pca$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord[, 1:4])

# Plot the correlation circle
a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = "PC1",  ylab = "PC2")
abline(h = 0, v = 0, lty = 2)
# Add active variables
arrows(0, 0, var.coord[, 1], var.coord[, 2], 
       length = 0.1, angle = 15, code = 2)
# Add labels
text(var.coord, labels=rownames(var.coord), cex = 1, adj=1)

#################################################################
#Make binary dataset based on residuals <=-1

resid_all_binary <- ifelse(resid_all<=-1, 1, 0)

#proportion of methylated cutsites
prop_methyl <- colSums(resid_all_binary) / nrow(resid_all_binary)
barplot(prop_methyl)

#MDS
resid_t_binary <- t(resid_all_binary)
d <- dist(resid_t_binary) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",  type="n")
text(x, y, labels = row.names(resid_t_binary), 
     col = rainbow(length(diam2))[rank(diam2)], cex=.7)


####################################################################
#DESeq2 analysis of differentially methylated loci
#NOTE: this does not seem to return useful information (it is skewed
#towards high count data, though maybe this could be dealt with)

library(DESeq2)
# Specify which columns are in which groups
sinfo <- read.delim("data3_sample_info.txt", header=FALSE)
deseq2.colData <- data.frame(condition=factor(sinfo[2:27,3]), 
                             type=factor(rep("paired-end", 26)))
rownames(deseq2.colData) <- colnames(data5)
deseq2.dds <- DESeqDataSetFromMatrix(countData = data5,
                                     colData = deseq2.colData, 
                                     design = ~ condition)
# Run Analysis
deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]

head(deseq2.res)

# Count number of hits with adjusted p-value less then 0.05
dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])

tmp <- deseq2.res
# The main plot
plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="darkgray",
     main="Diff. methylyated loci  (pval <= 0.05)",
     xlab="mean of normalized counts",
     ylab="Log2 Fold Change")

tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="chartreuse")
# 2 FC lines
abline(h=c(-1,1), col="blue")

#Get normalized counts and plot PCA
rld <- rlog(deseq2.dds)
plotPCA(rld, intgroup=c("condition", "type"))

counts2_cpm <- counts(deseq2.dds, normalized = TRUE)

sig_counts <- deseq_normcounts[c("11617", "445", "6418","6943", "791", "958"), ]

