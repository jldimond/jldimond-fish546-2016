## This script reads in a .str file from the ipyrad output and subsets the
## data to get only ddRAD data and extract only SNPs without missing data.

# Read in data file
data2.2 <- read.delim("data3.str", header=FALSE)
data <- data2.2[,colSums(is.na(data2.2))<nrow(data2.2)]
data2 <- t(data)
data3 <- as.data.frame(data2)
#Extract only ddRAD data 
data4 <- data3[,c(5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,47,48,51,52,55,56,59,60,63,64,69,70,73,74,77,78,81,82,85,86,89,90,93,94,97,98,101,102,105,106,109,110)]
#Remove any SNPs with missing data (-9 is the NA value)
data5 <- data4[!rowSums(data4[-c(1:2)] == -9) >= 1,]
data6 <-t(data5)
write.table(data6, file = "data3-2.str", row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################################
## This script reads in a .geno file from the ipyrad output and extracts ddRAD
## data without missing values.

# Read in .geno data file
a1 <- read.table("data3.geno", colClasses = 'character', header = FALSE)
a2 <- read.fwf("data3.geno", widths=rep(1, max(nchar(a1$V1))), colClasses = 'numeric', header=FALSE)
header <- read.delim("header_data3.txt", header=FALSE)
names <- t(header)
names2 <-as.vector(names)
colnames(a2) <- names2
#Select samples of interest (some have very low sample sizes)
a3 <- a2[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56)]
#Matrix with only ddr loci ***if including sample 101
a4 <- a3[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52)]
#Matrix with only ddr loci ***if excluding sample 101
a4 <- a3[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51)]
#Get rid of rows with any NAs (9)
a5 <- a4[!rowSums(a4[-c(1:2)] == 9) >= 1,]

#Read in sample info
sinfo <- read.table("sample_info.txt", colClasses = 'character', header = FALSE)
tsinfo <- t(sinfo)
type <- tsinfo[3,]
diam <- tsinfo[5,]

#MDS
a6 <- t(a5)
d <- dist(a6) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS")
text(x, y, labels = row.names(a6), col = rainbow(length(diam))[rank(diam)], cex=.8)

sinfo <- read.table("sample_info.txt", colClasses = 'character', header = FALSE)
tsinfo <- t(sinfo)
type <- tsinfo[3,]
