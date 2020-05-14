####################################################################################################################
# take RNA-seq data from:
# Desai, N. et al. (2020). 	Spectrum of Viral Load and Host Response Seen in Autopsies of SARS-CoV-2 Infected Lungs.
# Retrieved on May 14, 2020 from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150316
# find differentially expressed genes in COVID-19 patients
# make a heatmap of significant differentially expressed genes
####################################################################################################################

library(DESeq2)
library(limma)
library(gplots)
library(RColorBrewer)
library(edgeR)
library(DEFormats)

#define functions for coloring and distance
redgreen <- function(n) {
  c(hsv(h=2/6, v=seq(1,0,length=n/2) ), hsv(h=0/6, v=seq(0,1,length=n/2)) )
}
corr.dist=function(x) { as.dist(1-cor(t(x))) }

#extract count data and get counts per million
countdata <- read.table(file="GSE150316_RawCounts_Final.txt")
cpm <- cpm(countdata)
#establish a threshold
thresh <- cpm > 0.5
#only keep counts above threshold and 
#with more than one sample over that threshold
keep <- rowSums(thresh) >= 2
count.thresh <- countdata[keep,]
#add a disease identifier
names(count.thresh) <- names(count.thresh)<-
  paste(ifelse(substr(names(count.thresh),1,3) =="Neg","H","D"),names(count.thresh),sep = ":")
#make a DGEList object
dgeObj <- DGEList(count.thresh)
#take the regularized log of the counts in the DGE object
rld <- rlog(as.matrix(dgeObj$counts))
#set the gene names for the regularized log counts 
#factor samples by whether or not they have the disease
rownames(rld) <- rownames(dgeObj$counts)
dgeObj[["samples"]][["group"]] <- c(rep(1,32),rep(0,5))
dgeObj[["samples"]][["group"]] <- factor(dgeObj[["samples"]][["group"]],
                                         unique(dgeObj[["samples"]][["group"]]))
#convert to dds object
dds <- as.DESeqDataSet(dgeObj)

#get differential expression results
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
#select genes with a log fold change over 1 
#and an adjusted pvalue below 0.05
slct <- res@rownames[abs(res$log2FoldChange)>1 & res$padj < 0.05]
slct.genes <- slct[!is.na(slct)]

#define major clusters of genes
#color each cluster
hc.de <- hclust(corr.dist(rld[slct.genes,]),method = "ward.D2")
clusters.de <- cutree(hc.de,h=2)
rcol.de = brewer.pal(3,"Paired")[clusters.de]

#make a heatmap based on these genes
heatmap.2(rld[slct.genes,],trace="none",
          labRow=F,col=redgreen(100),scale="row",
          margins=c(7,7), distfun = corr.dist,
          ColSideColors = c(rep("red",32),rep("green",5)),
          RowSideColors=rcol.de,
          hclustfun=function(x) { hclust(x,method="ward.D2")},
          main="SARS-Cov-2 Gene Clusters")
legend(.001,.8,legend=1:3,pch=15,col=brewer.pal(3,"Paired"))

