# load/install limma, Glimma, edgeR, Mus.musculus
# then
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)


## Data packaging

# Reading in count-data

# download data (.tar) and extract files
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="data/GSE63310_RAW.tar", mode="wb")
utils::untar("data/GSE63310_RAW.tar", exdir = "./data")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")

# to read raw gene-level counts for the given 9 samples
files <- c("data/GSM1545535_10_6_5_11.txt", "data/GSM1545536_9_6_5_11.txt", 
   "data/GSM1545538_purep53.txt", "data/GSM1545539_JMS8-2.txt", 
   "data/GSM1545540_JMS8-3.txt", "data/GSM1545541_JMS8-4.txt", 
   "data/GSM1545542_JMS8-5.txt", "data/GSM1545544_JMS9-P7c.txt", 
   "data/GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)		# the 1st file is tab-separated and of only one line

# read the text files of raw gene-level counts and combine them them into a table of counts
# The resulting DGEList-object contains the table of counts, each row associated with unique Entrez gene identifiers (IDs) and
# nine columns associated with the individual samples
x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)	# returns the dimensions of the object x (table)

# Organize sample information
# each of the 9 text files contain a unique sample, which is of one of the 3 cell types (basal, LP and ML)
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples


## Data pre-processing

# Transformations from the raw counts to counts per million (CPM) and log2-counts per million (log-CPM)
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# exclude the genes that are not expressed (zero counts)
table(rowSums(x$counts==0)==9)

# CPM = 1 is used as the expression cutoff (used to separate expressed genes from unexpressed genes)
# CPM = 1 means that a gene is expressed if it has at least 20 counts in the sample with the lowest sequencing depth (JMS9-P8c, library size approx. 20 million) or 
# at least 76 counts in the sample with the greatest sequencing depth (JMS8-3, library size approx. 76 million). 
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# plot the cutoff data
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
 den <- density(lcpm[,i])
 lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
   den <- density(lcpm[,i])
   lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# Normalize gene expression distributions

# norm factors of each sample (trimmed mean of M-values (TMM))
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# visualize the distributions
# to give good visual representation
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

# distributions of pre-normalized data
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

# distributions of normalized data
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

# Unsupervised clustering of samples

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")


## Differential expression analysis

# Creating a design matrix and contrasts
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

# Remove heteroscedascity from count data
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

# Examine the number of DE genes
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

# Examine individual DE genes
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)

# graphical representations of differential expression results
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))


