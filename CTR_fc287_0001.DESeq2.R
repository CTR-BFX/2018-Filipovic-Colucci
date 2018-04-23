#!/usr/local/bin/Rscript

#
# CTR_fc287_0001 ::: HTSeq to DESeq2 Analysis
#
# 
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#

#
# initial install of packages
#
#source("http://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
#devtools::install_github('slowkow/ggrepel')
#install.packages("ggrepel")
#install.packages("gplots")

library('DESeq2')
library('ggplot2')
library('RColorBrewer')
library("cowplot")
library("pheatmap")
library("ggrepel")
library("reshape2")
library("biomaRt")
library("matrixStats")
library("plyr")
library("UpSetR")
library("genefilter")
library("scales")
library("gtools")
library("Cairo")

message("+-------------------------------------------------------------------------------")
message("+ Set up some constants e.g. base directories")
message("+-------------------------------------------------------------------------------")

Project  <- "CTR_fc287_0001"
Base.dir <- "/Users/rhamilto/Documents/CTR-Groups/Francesco_Colucci/CTR_fc287_0001"
setwd(Base.dir)
HTSeq.dir <- paste(Base.dir,"/HTSeq_Counts", sep="")

elementTextSize <- 10

message("+-------------------------------------------------------------------------------")
message("+ Set up the sample table")
message("+-------------------------------------------------------------------------------")
#1   u-cNK-1     L2V18DRBC01
#2   u-trNK-1    L2V18DRBC02
#3   u-ILC1-1    L2V18DRBC03
#4   l-ILC1-1    L2V18DRBC04
#5   l-cNK-1     L2V18DRBC05

#6   u-cNK-2     L2V18DRBC06
#7   u-trNK-2    L2V18DRBC07
#8   u-ILC1-2    L2V18DRBC08
#9   l-ILC1-2    L2V18DRBC09
#10  l-cNK-2     L2V18DRBC10

#11  u-cNK-3     L2V18DRBC11
#12  u-trNK-3    L2V18DRBC12
#13  u-ILC1-3    L2V18DRBC13
#14  l-ILC1-3    L2V18DRBC14
#15  l-cNK-3     L2V18DRBC15

#u=uterus
#l=liver
#c=conventional
#tr=tissue resident
#ILC1=innate lymphoid cell
#NK=natural killer

sampleFiles      <- grep('*counts.txt',list.files(HTSeq.dir),value=TRUE)
sampleNames      <- gsub(".HGHVHBBXX.s_1.r_1_trimmed_tophat.accepted_hits_merged.bam_htseq_counts.txt", "", sampleFiles)
sampleBarcodes   <- gsub("SLX-9343.", "", sampleNames)
sampleBarcodes   <- gsub("L2V18DRBC", "", sampleBarcodes)

sampleCell       <- c("cNK", "trNK", "ILC1", "ILC1", "cNK",                "cNK", "trNK", "ILC1", "ILC1", "cNK",               "cNK", "trNK", "ILC1", "ILC1", "cNK")
sampleTissue     <- c("uterus","uterus","uterus","liver","liver",          "uterus","uterus","uterus","liver","liver",         "uterus","uterus","uterus","liver","liver")
sampleCollated1  <- c("none","utrNK_uILC1","utrNK_uILC1","lILC1", "none",  "none","utrNK_uILC1","utrNK_uILC1","lILC1","none",  "none","utrNK_uILC1","utrNK_uILC1", "lILC1","none")
sampleCollated2  <- c("ucNK","uILC1_utrNK","uILC1_utrNK","none","none",    "ucNK","uILC1_utrNK","uILC1_utrNK","none","none",   "ucNK","uILC1_utrNK","uILC1_utrNK","none","none")

sampleGroup      <- factor(paste0(sampleCell, sampleTissue))

sampleTable      <- data.frame(sampleName=sampleBarcodes, fileName=sampleFiles, cell=sampleCell, tissue=sampleTissue, group=sampleGroup, collated1=sampleCollated1, collated2=sampleCollated2)
print(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Removing samples identified by QC")
message("+-------------------------------------------------------------------------------")
sampleTable <- sampleTable[-c(3,11),]
print(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Retrieve ensEMBL annotations")
message("+-------------------------------------------------------------------------------")
ensembl    <- useEnsembl(biomart="ensembl", host="http://mar2016.archive.ensembl.org", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene', 'chromosome_name', 'gene_biotype'), mart = ensembl)  
head(ensEMBL2id)
nrow(ensEMBL2id)

write.csv2(ensEMBL2id, file=paste0(Project, "_ensEMBL2id.csv"))
ensEMBL2id       <- read.table(paste0(Project, "_ensEMBL2id.csv"), sep=";", header=TRUE, stringsAsFactors=FALSE)
head(ensEMBL2id)
nrow(ensEMBL2id)


message("+-------------------------------------------------------------------------------")
message("+ Retrieve average transcript lengths")
message("+-------------------------------------------------------------------------------")

# From: http://seqanswers.com/forums/archive/index.php/t-39797.html

library(GenomicRanges)
library(rtracklayer)

GTFfile                               <- "/usr/local/Genomes/Mus_musculus/GRCm38/Mus_musculus.GRCm38.84.gtf"
GTF                                   <- import.gff(GTFfile, format="gtf", genome="GRCm38.84", feature.type="exon")
grl                                   <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF                            <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id   <- rep(names(grl), elementNROWS(grl))
elementMetadata(reducedGTF)$widths    <- width(reducedGTF)
calc_length                           <- function(x) { sum(elementMetadata(x)$widths) }
ensembl.aveTranscLens                 <- as.data.frame(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
colnames(ensembl.aveTranscLens)       <- c("Length")
ensembl.aveTranscLens$ensembl_gene_id <- rownames(ensembl.aveTranscLens)

ensEMBL2id                            <- merge(ensEMBL2id, ensembl.aveTranscLens, by="ensembl_gene_id" )
head(ensEMBL2id)
nrow(ensEMBL2id)



message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")
ddsHTSeq.group      <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ group) 
ddsHTSeq.tissuecell <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ tissue + cell)
ddsHTSeq.collated1  <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ collated1)
ddsHTSeq.collated2  <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ collated2)

ddsHTSeq.tissue     <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ tissue)




message("+-------------------------------------------------------------------------------")
message("+ Create dds objects")
message("+-------------------------------------------------------------------------------")

dds.group      <- DESeq(ddsHTSeq.group)
dds.tissuecell <- DESeq(ddsHTSeq.tissuecell)
dds.collated1  <- DESeq(ddsHTSeq.collated1)
dds.collated2  <- DESeq(ddsHTSeq.collated2)
dds.tissue     <- DESeq(ddsHTSeq.tissue)


message("+-------------------------------------------------------------------------------")
message("+ Calculate FPKMs and TPMs from dds")
message("+-------------------------------------------------------------------------------")

# Remove duplicates from the annotation table
ensEMBL2id <- ensEMBL2id[ !duplicated(ensEMBL2id$ensembl_gene_id), ]
head(ensEMBL2id)
nrow(ensEMBL2id)

# make sure they are in the same order as the dds object
ensEMBL2id$ensembl_gene_id <- ordered(ensEMBL2id$ensembl_gene_id, rownames(ddsHTSeq.group) ) 
head(ensEMBL2id)
nrow(ensEMBL2id)

mcols(dds.group)$basepairs <- ensEMBL2id$Length
dds.group.FPKM <- fpkm(dds.group)
head(dds.group.FPKM)
write.csv(dds.group.FPKM, file=paste(Project, "_DESeq2NormCounts_2_FPKM.csv", sep=""))


#TPMi=( FPKMi / sum(FPKMj ) * 10^6
dds.group.TPM <- do.call(cbind, lapply(1:ncol(dds.group.FPKM), function(i) { (dds.group.FPKM[,i] / sum(dds.group.FPKM[,i])) * 1e6}))
colnames(dds.group.TPM) <- sampleTable$sampleName
head(dds.group.TPM)
write.csv(dds.group.TPM, file=paste(Project, "_DESeq2NormCounts_2_FPKM_2TPM.csv", sep=""))



message("+-------------------------------------------------------------------------------")
message("+ Create results object")
message("+-------------------------------------------------------------------------------")

normCounts.group                     <- counts(dds.group,      normalized=TRUE)
normCounts.tissuecell                <- counts(dds.tissuecell, normalized=TRUE)
normCounts.collated1                 <- counts(dds.collated1,  normalized=TRUE)
normCounts.collated2                 <- counts(dds.collated2,  normalized=TRUE)


normCounts.group.df                  <- as.data.frame(normCounts.group)
normCounts.group.df$ensembl_gene_id  <- rownames(normCounts.group.df)
normCounts.group.df.annot            <- merge(normCounts.group.df,ensEMBL2id, by="ensembl_gene_id" )
head(normCounts.group.df.annot)

write.csv(normCounts.group.df.annot, file=paste(Project, "NormalisedCounts_normCounts.group.csv", sep=""))


significance=0.05
foldchange=2


message("+-------------------------------------------------------------------------------")
message("+ Perform Comparisons of defined groups, tissues, cells")
message("+-------------------------------------------------------------------------------")

res_uterus_vs_liver                        <- results(dds.tissuecell, contrast=c("tissue", "uterus", "liver"))

res_uterus_cNK_vs_uterus_trNK              <- results(dds.group,      contrast=c("group", "cNKuterus",  "trNKuterus")) # 7a
res_uterus_trNK_vs_uterus_cNK              <- results(dds.group,      contrast=c("group", "trNKuterus", "cNKuterus" )) # 7b

res_uterus_cNK_vs_uterus_ILC1              <- results(dds.group,      contrast=c("group", "cNKuterus",  "ILC1uterus")) # 6
res_uterus_ILC1_vs_uterus_trNK             <- results(dds.group,      contrast=c("group", "ILC1uterus", "trNKuterus")) # 10

res_liver_cNK_vs_liver_ILC1                <- results(dds.group,      contrast=c("group", "cNKliver",   "ILC1liver" )) # 2

res_uterus_cNK_vs_liver_cNK                <- results(dds.group,      contrast=c("group", "cNKuterus",  "cNKliver"  )) # 1
res_uterus_ILC1_vs_liver_ILC1              <- results(dds.group,      contrast=c("group", "ILC1uterus", "ILC1liver" )) # 8

res_uterus_trNK_vs_liver_ILC1              <- results(dds.group,      contrast=c("group", "trNKuterus", "ILC1liver" )) # 9
res_uterus_trNK_vs_liver_cNK               <- results(dds.group,      contrast=c("group", "trNKuterus", "cNKliver"  )) # 4

res_liver_cNK_vs_liver_ILC1                <- results(dds.group,      contrast=c("group", "cNKliver",   "ILC1liver" )) # 5
res_uterus_cNK_vs_liver_ILC1               <- results(dds.group,      contrast=c("group", "cNKuterus",  "ILC1liver" )) # 3

res_uterus_trNK_uterus_ILC1_vs_liver_ILC1  <- results(dds.collated1,  contrast=c("collated1", "utrNK_uILC1", "lILC1"))
res_uterus_ILC1_uterus_trNK_vs_uterus_cNK  <- results(dds.collated2,  contrast=c("collated2", "uILC1_utrNK", "ucNK"))

message("+-------------------------------------------------------------------------------")
message("+ Filter and output comparisons of defined groups, tissues, cells")
message("+-------------------------------------------------------------------------------")


functionGetSigDEG <- function(result, Project, Title, significance, foldchange, ensEMBL2id, Counts) {
  # Filter by adjusted p-value
  result.df <- as.data.frame(subset(result, padj <= significance  & abs(log2FoldChange) >= foldchange))  
  # Print out the number of significant hits
  message(paste("+   ", Title, "            = ", nrow(result.df),sep=""))
  # Annotate Genes
  result.df$ensembl_gene_id <- rownames(result.df)
  # Annotate from the live ensEMBL table generated above
  result.ann <- merge(result.df,ensEMBL2id, by="ensembl_gene_id")
  # Tidy up the description field
  result.ann$description            <- gsub("..Source.*", "", result.ann$description)
  # Write the results to file
  write.csv(result.ann[order(abs(result.ann$log2FoldChange),decreasing=TRUE),], file=paste(Project, "_DESeq2_DEGs_sig", significance, '_fc', foldchange, '_', Title, '.ann.csv', sep=""))
  
  result.ann.rlog <- merge(result.ann, Counts, by="ensembl_gene_id")
  
  result.ann.rlog <- result.ann.rlog[,c(2:25)]
  
  write.csv(result.ann.rlog[order(abs(result.ann.rlog$log2FoldChange),decreasing=TRUE),], file=paste(Project, "_DESeq2_DEGs_sig", significance, '_fc', foldchange, '_', Title, '.ann.rlog.csv', sep=""))
  
  return(result.ann)
  #return(result.ann.rlog)
}


rld.group                       <- rlogTransformation(dds.group,      blind=T)
rld.group.assay                 <- as.data.frame( assay(rld.group) )
rld.group.assay$ensembl_gene_id <- rownames(rld.group.assay)
rld.group.ann                   <- merge(rld.group.assay, ensEMBL2id, by="ensembl_gene_id")

res_uterus_vs_liver.ann                        <- functionGetSigDEG(as.data.frame(res_uterus_vs_liver),                       Project,"res_uterus_vs_liver",                       significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))

res_uterus_cNK_vs_uterus_trNK.ann              <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_uterus_trNK),             Project,"res_uterus_cNK_vs_uterus_trNK",             significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))

res_uterus_trNK_vs_uterus_cNK.ann               <- functionGetSigDEG(as.data.frame(res_uterus_trNK_vs_uterus_cNK),             Project,"res_uterus_trNK_vs_uterus_cNK",             significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))

res_uterus_cNK_vs_uterus_ILC1.ann              <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_uterus_ILC1),             Project,"res_uterus_cNK_vs_uterus_ILC1",             significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_ILC1_vs_uterus_trNK.ann             <- functionGetSigDEG(as.data.frame(res_uterus_ILC1_vs_uterus_trNK),            Project,"res_uterus_ILC1_vs_uterus_trNK",            significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_liver_cNK_vs_liver_ILC1.ann                <- functionGetSigDEG(as.data.frame(res_liver_cNK_vs_liver_ILC1),               Project,"res_liver_cNK_vs_liver_ILC1",               significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_cNK_vs_liver_cNK.ann                <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_liver_cNK),               Project,"res_uterus_cNK_vs_liver_cNK",               significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_ILC1_vs_liver_ILC1.ann              <- functionGetSigDEG(as.data.frame(res_uterus_ILC1_vs_liver_ILC1),             Project,"res_uterus_ILC1_vs_liver_ILC1",             significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_trNK_vs_liver_ILC1.ann              <- functionGetSigDEG(as.data.frame(res_uterus_trNK_vs_liver_ILC1),             Project,"res_uterus_trNK_vs_liver_ILC1",             significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_trNK_vs_liver_cNK.ann               <- functionGetSigDEG(as.data.frame(res_uterus_trNK_vs_liver_cNK),              Project,"res_uterus_trNK_vs_liver_cNK",              significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_liver_cNK_vs_liver_ILC1.ann                <- functionGetSigDEG(as.data.frame(res_liver_cNK_vs_liver_ILC1),               Project,"res_liver_cNK_vs_liver_ILC1",               significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_cNK_vs_liver_ILC1.ann               <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_liver_ILC1),              Project,"res_uterus_cNK_vs_liver_ILC1",              significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))

res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.ann  <- functionGetSigDEG(as.data.frame(res_uterus_trNK_uterus_ILC1_vs_liver_ILC1), Project,"res_uterus_trNK_uterus_ILC1_vs_liver_ILC1", significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.ann  <- functionGetSigDEG(as.data.frame(res_uterus_ILC1_uterus_trNK_vs_uterus_cNK), Project,"res_uterus_ILC1_uterus_trNK_vs_uterus_cNK", significance,foldchange,ensEMBL2id, as.data.frame(rld.group.ann))





#head(res_uterus_vs_liver.ann[order(abs(res_uterus_vs_liver.ann$log2FoldChange),decreasing=TRUE),], 20)
#test <- as.data.frame( head(results(dds.tissuecell, contrast=c("tissue", "uterus", "liver")), 20) )
#for(i in rownames(test)){ print(i) }

message("+-------------------------------------------------------------------------------")
message("+ Create MA Plots")
message("+-------------------------------------------------------------------------------")



favoritegenesA <- NULL
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000021477") #Ctsl
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000038642") #Ctss
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000021939") #Ctsb
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000024621") #Csf1r
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000073421") #H2-Ab1
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000036594") #H2-Aa
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000024610") #Cd74
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000046805") #Mpeg1
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000069516") #Lyz2
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000004207") #Psap
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000019122") #Ccl9

favoritegenesB <- NULL
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000023132") #Gzma	
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000015441") #Gzmf
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000059256") #Gzmd
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000040284") #Gzmg
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000022156") #Gzme
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000079186") #Gzmc

favoritegenesC <- NULL
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000029304") #Spp1
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000035493") #Tgfbi
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000040552") #C3ar1

favoritegenesD <- NULL
favoritegenesD <- append(favoritegenesD, "ENSMUSG00000035042") #Ccl5


functionCustomMAPlot <- function(results, Project, Title, significance, log2FoldChange, favoritegenesA, favouritegenesB, favouritegenesC, favouritegenesD) {
  # Get annotation infor for a favorite set of ensEMBL ids
  labeldata.rows            <- match(favoritegenesA, row.names(results))
  labeldata                 <- results[labeldata.rows,]
  labeldata$ensembl_gene_id <- rownames(labeldata)
  labeldata.ann             <- merge(labeldata, ensEMBL2id, by="ensembl_gene_id")
  
  labeldataB.rows            <- match(favoritegenesB, row.names(results))
  labeldataB                 <- results[labeldataB.rows,]
  labeldataB$ensembl_gene_id <- rownames(labeldataB)
  labeldataB.ann             <- merge(labeldataB, ensEMBL2id, by="ensembl_gene_id")
  
  labeldataC.rows            <- match(favoritegenesC, row.names(results))
  labeldataC                 <- results[labeldataC.rows,]
  labeldataC$ensembl_gene_id <- rownames(labeldataC)
  labeldataC.ann             <- merge(labeldataC, ensEMBL2id, by="ensembl_gene_id")
  
  labeldataD.rows            <- match(favoritegenesD, row.names(results))
  labeldataD                 <- results[labeldataD.rows,]
  labeldataD$ensembl_gene_id <- rownames(labeldataD)
  labeldataD.ann             <- merge(labeldataD, ensEMBL2id, by="ensembl_gene_id")
  
  results$log2FoldChange[results$log2FoldChange > 12.5]  <- 12.5
  results$log2FoldChange[results$log2FoldChange < -12.5] <- -12.5
  
  plt <- ggplot(data = results, aes(x=baseMean, y=log2FoldChange )) + 
    geom_abline(intercept = log2FoldChange, slope = 0, colour='red', alpha=0.25) + 
    geom_abline(intercept = -log2FoldChange, slope = 0, colour='red', alpha=0.25) +
    geom_point(size=0.5, alpha=0.5, col="black") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange >= 0)), size=1, alpha=0.5,  col="red") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange < 0)),  size=1, alpha=0.5,  col="blue") +
    
    geom_point(data=subset(labeldata.ann, log2FoldChange > 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkred',  shape=21, stroke=0.5) +
    geom_point(data=subset(labeldata.ann, log2FoldChange < 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkblue', shape=21, stroke=0.5) +
    
    geom_point(data=labeldataB.ann, aes( x=baseMean, y=log2FoldChange), size=3, alpha=1.0, color='pink', shape=21, stroke=0.5) +
    geom_point(data=labeldataB.ann, aes( x=baseMean, y=log2FoldChange), size=1, alpha=1.0, color='pink', stroke=0.5) +
    
    geom_point(data=labeldataC.ann, aes( x=baseMean, y=log2FoldChange), size=3, alpha=1.0, color='orange', shape=21, stroke=0.5) +
    geom_point(data=labeldataC.ann, aes( x=baseMean, y=log2FoldChange), size=1, alpha=1.0, color='orange', stroke=0.5) +
    
    geom_point(data=labeldataD.ann, aes( x=baseMean, y=log2FoldChange), size=3, alpha=1.0, color='blue', shape=21, stroke=0.5) +
    geom_point(data=labeldataD.ann, aes( x=baseMean, y=log2FoldChange), size=1, alpha=1.0, color='blue', stroke=0.5) +
    
    
    geom_label_repel(data=subset(labeldata.ann, log2FoldChange > 0), 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='red', colour='white', point.padding = unit(0.25, "lines"),  size=6, segment.color = 'darkred',  nudge_x = 1) +
    geom_label_repel(data=subset(labeldata.ann, log2FoldChange < 0), 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='blue', colour='white', point.padding = unit(0.25, "lines"), size=6, segment.color = 'darkblue', nudge_x = -1) +
    
    geom_label_repel(data=labeldataB.ann, 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='pink', colour='white', point.padding = unit(0.25, "lines"),  size=6, segment.color = 'pink',  nudge_x = -2) +
    
    geom_label_repel(data=labeldataC.ann, 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='orange', colour='white', point.padding = unit(0.25, "lines"),  size=6, segment.color = 'orange',  nudge_x = -3) +
    
    geom_label_repel(data=labeldataD.ann, 
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name), 
                     fill='blue', colour='white', point.padding = unit(0.25, "lines"),  size=6, segment.color = 'blue',  nudge_x = 2) +
    
    scale_x_log10() + 
    #scale_y_reverse() +
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") + ggtitle(paste(Project, " DESeq2 MA ", Title, " [fc ", log2FoldChange, ", sig ", significance, "]", sep="")) 
  
  pdf(paste(Project, "_DESeq_MA_", Title, "_fc", log2FoldChange, "_sig", significance, ".pdf", sep=""),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt })
  dev.off()
}


#functionCustomMAPlot(as.data.frame(res_uterus_vs_liver),                         Project, "res_uterus_vs_liver",                        significance, foldchange, favoritegenes)
functionCustomMAPlot(as.data.frame(res_uterus_cNK_vs_uterus_trNK),               Project, "res_uterus_cNK_vs_uterus_trNK",              significance, foldchange, favoritegenesA, favoritegenesB, favouritegenesC, favouritegenesD)
functionCustomMAPlot(as.data.frame(res_uterus_trNK_vs_uterus_cNK),               Project, "res_uterus_trNK_vs_uterus_cNK",              significance, foldchange, favoritegenesA, favoritegenesB, favouritegenesC, favouritegenesD)

#functionCustomMAPlot(as.data.frame(res_uterus_cNK_vs_uterus_ILC1),               Project, "res_uterus_cNK_vs_uterus_ILC1",              significance, foldchange, favoritegenes)
#functionCustomMAPlot(as.data.frame(res_uterus_ILC1_vs_uterus_trNK),              Project, "res_uterus_ILC1_vs_uterus_trNK",             significance, foldchange, favoritegenes)
#functionCustomMAPlot(as.data.frame(res_liver_cNK_vs_liver_ILC1),                 Project, "res_liver_cNK_vs_liver_ILC1",                significance, foldchange, favoritegenes)
#functionCustomMAPlot(as.data.frame(res_uterus_cNK_vs_liver_cNK),                 Project, "res_uterus_cNK_vs_liver_cNK",                significance, foldchange, favoritegenes)
#functionCustomMAPlot(as.data.frame(res_uterus_ILC1_vs_liver_ILC1),               Project, "res_uterus_ILC1_vs_liver_ILC1",              significance, foldchange, favoritegenes)

#functionCustomMAPlot(as.data.frame(res_uterus_trNK_vs_liver_ILC1),               Project, "res_uterus_trNK_vs_liver_ILC1",              significance, foldchange, favoritegenes)
#functionCustomMAPlot(as.data.frame(res_uterus_trNK_vs_liver_cNK),                Project, "res_uterus_trNK_vs_liver_cNK",               significance, foldchange, favoritegenes)

#functionCustomMAPlot(as.data.frame(res_uterus_trNK_uterus_ILC1_vs_liver_ILC1),   Project, "res_uterus_trNK_uterus_ILC1_vs_liver_ILC1",  significance, foldchange, favoritegenes)
#functionCustomMAPlot(as.data.frame(res_uterus_ILC1_uterus_trNK_vs_uterus_cNK),   Project, "res_uterus_ILC1_uterus_trNK_vs_uterus_cNK",  significance, foldchange, favoritegenes)



message("+-------------------------------------------------------------------------------")
message("+ UpSetR")
message("+-------------------------------------------------------------------------------")

listInput <- list( 'uterus vs liver'                        = res_uterus_vs_liver.ann$ensembl_gene_id, 
                   'uterus cNK vs uterus trNK'              = res_uterus_cNK_vs_uterus_trNK.ann$ensembl_gene_id, 
                   'uterus cNK vs uterus ILC1'              = res_uterus_cNK_vs_uterus_ILC1.ann$ensembl_gene_id,
                   'uterus ILC1 vs uterus trNK'             = res_uterus_ILC1_vs_uterus_trNK.ann$ensembl_gene_id,
                   'liver cNK vs liver ILC1'                = res_liver_cNK_vs_liver_ILC1.ann$ensembl_gene_id,
                   'uterus cNK vs liver cNK'                = res_uterus_cNK_vs_liver_cNK.ann$ensembl_gene_id,
                   'uterus ILC1 vs liver ILC1'              = res_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                   'uterus trNK vs liver ILC1'              = res_uterus_trNK_vs_liver_ILC1.ann$ensembl_gene_id,
                   'uterus trNK vs liver cNK'               = res_uterus_trNK_vs_liver_cNK.ann$ensembl_gene_id,
                   'uterus trNK and uterus ILC1 vs liver ILC1'  = res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                   'uterus ILC1 and uterus trNK vs uterus cNK'  = res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.ann$ensembl_gene_id
                 )


listInput <- list("uterus vs liver"                         = res_uterus_vs_liver.ann$ensembl_gene_id,
                  'uterine trNK vs uterine cNK'             = res_uterus_cNK_vs_uterus_trNK.ann$ensembl_gene_id,
                  'uterine ILC1 vs uterine cNK'             = res_uterus_cNK_vs_uterus_ILC1.ann$ensembl_gene_id,
                  'uterine ILC1 vs uterine trNK'            = res_uterus_ILC1_vs_uterus_trNK.ann$ensembl_gene_id,
                  'liver ILC1 vs liver cNK'                 = res_liver_cNK_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine cNK vs liver cNK'                = res_uterus_cNK_vs_liver_cNK.ann$ensembl_gene_id,
                  'uterine ILC1 vs liver ILC1'              = res_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine trNK vs liver ILC1'              = res_uterus_trNK_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine trNK vs liver cNK'               = res_uterus_trNK_vs_liver_cNK.ann$ensembl_gene_id,
                  'uterine trNK and ILC1 vs liver ILC1'     = res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine trNK and ILC1 vs uterine cNK'    = res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.ann$ensembl_gene_id
                 )

pdf(paste(Project, "_DESeq_UpSetR.pdf", sep=""),width=10,height=7, onefile=FALSE)
par(bg=NA)
upset(fromList(listInput), nsets = 11, sets = rev(colnames(fromList(listInput))), keep.order = TRUE, empty.intersections = "on", 
      nintersects = 11, sets.bar.color = rev(c("blue", "blue", "blue", "green", "green", "red", "red", "red", "red", "red", "red")), 
      sets.x.label = "Number of differentially expressed genes", mainbar.y.label = "Intersections of Differentially Expressed Genes")
dev.off()



L01 <- res_uterus_vs_liver.ann$ensembl_gene_id
L02 <- res_uterus_cNK_vs_uterus_trNK.ann$ensembl_gene_id 
L03 <- res_uterus_cNK_vs_uterus_ILC1.ann$ensembl_gene_id
L04 <- res_uterus_ILC1_vs_uterus_trNK.ann$ensembl_gene_id
L05 <- res_liver_cNK_vs_liver_ILC1.ann$ensembl_gene_id
L06 <- res_uterus_cNK_vs_liver_cNK.ann$ensembl_gene_id
L07 <- res_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id
L08 <- res_uterus_trNK_vs_liver_ILC1.ann$ensembl_gene_id
L09 <- res_uterus_trNK_vs_liver_cNK.ann$ensembl_gene_id
L10 <- res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id
L11 <- res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.ann$ensembl_gene_id


L01.uniq           <- as.data.frame(setdiff( L01, c(L02,L03,L04,L05,L06,L07,L08,L09,L10,L11) ) ) 
nrow(L01.uniq)
colnames(L01.uniq) <- "ensembl_gene_id"
L01.uniq.ann <- merge(L01.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L01.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_vs_liver.ann.csv", sep=""))

L02.uniq           <- as.data.frame(setdiff( L02, c(L01,L03,L04,L05,L06,L07,L08,L09,L10,L11) ) ) 
nrow(L02.uniq)
colnames(L02.uniq) <- "ensembl_gene_id"
L02.uniq.ann <- merge(L02.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L02.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_cNK_vs_uterus_trNK.ann.csv", sep=""))

L03.uniq           <- as.data.frame(setdiff( L03, c(L01,L02,L04,L05,L06,L07,L08,L09,L10,L11) ) ) 
nrow(L03.uniq)
colnames(L03.uniq) <- "ensembl_gene_id"
L03.uniq.ann <- merge(L03.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L03.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_cNK_vs_uterus_ILC1.ann.csv", sep=""))

L04.uniq           <- as.data.frame(setdiff( L04, c(L01,L02,L03,L05,L06,L07,L08,L09,L10,L11) ) ) 
nrow(L04.uniq)
colnames(L04.uniq) <- "ensembl_gene_id"
L04.uniq.ann <- merge(L04.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L04.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_ILC1_vs_uterus_trNK.ann.csv", sep=""))

L05.uniq           <- as.data.frame(setdiff( L05, c(L01,L02,L03,L04,L06,L07,L08,L09,L10,L11) ) ) 
nrow(L05.uniq)
colnames(L05.uniq) <- "ensembl_gene_id"
L05.uniq.ann <- merge(L05.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L05.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_liver_cNK_vs_liver_ILC1.ann.csv", sep=""))

L06.uniq           <- as.data.frame(setdiff( L06, c(L01,L02,L03,L04,L05,L07,L08,L09,L10,L11) ) ) 
nrow(L06.uniq)
colnames(L06.uniq) <- "ensembl_gene_id"
L06.uniq.ann <- merge(L06.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L06.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_cNK_vs_liver_cNK.ann.csv", sep=""))

L07.uniq           <- as.data.frame(setdiff( L07, c(L01,L02,L03,L04,L05,L06,L08,L09,L10,L11) ) ) 
nrow(L07.uniq)
colnames(L07.uniq) <- "ensembl_gene_id"
L07.uniq.ann <- merge(L07.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L07.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_ILC1_vs_liver_ILC1.ann.csv", sep=""))

L08.uniq           <- as.data.frame(setdiff( L08, c(L01,L02,L03,L04,L05,L06,L07,L09,L10,L11) ) ) 
nrow(L08.uniq)
colnames(L08.uniq) <- "ensembl_gene_id"
L08.uniq.ann <- merge(L08.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L08.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_trNK_vs_liver_ILC1.ann.csv", sep=""))

L09.uniq           <- as.data.frame(setdiff( L09, c(L01,L02,L03,L04,L05,L06,L07,L08,L10,L11) ) ) 
nrow(L09.uniq)
colnames(L09.uniq) <- "ensembl_gene_id"
L09.uniq.ann <- merge(L09.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L09.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_trNK_vs_liver_cNK.ann.csv", sep=""))

L10.uniq           <- as.data.frame(setdiff( L10, c(L01,L02,L03,L04,L05,L06,L07,L08,L09,L11) ) ) 
nrow(L10.uniq)
colnames(L10.uniq) <- "ensembl_gene_id"
L10.uniq.ann <- merge(L10.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L10.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.ann.csv", sep=""))

L11.uniq           <- as.data.frame(setdiff( L11, c(L01,L02,L03,L04,L05,L06,L07,L08,L09,L10) ) ) 
nrow(L11.uniq)
colnames(L11.uniq) <- "ensembl_gene_id"
L11.uniq.ann <- merge(L11.uniq, ensEMBL2id, by="ensembl_gene_id")
write.csv(L11.uniq.ann, file=paste(Project, "_UpSetR_Uniq_res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.ann.csv", sep=""))

message("+-------------------------------------------------------------------------------")
message("+ Run transformations (rlog transforms)")
message("+-------------------------------------------------------------------------------")

rld.group                       <- rlogTransformation(dds.group,      blind=T)
rld.group.assay                 <- as.data.frame( assay(rld.group) )
rld.group.assay$ensembl_gene_id <- rownames(rld.group.assay)
rld.group.ann                   <- merge(rld.group.assay, ensEMBL2id, by="ensembl_gene_id")

write.table(rld.group.ann, file=paste(Project, "_rlogTransformed.tsv", sep=""), sep="\t")

rld.group.ann.proteinCoding <- rld.group.ann[rld.group.ann$gene_biotype=='protein_coding',]
write.table(rld.group.ann.proteinCoding, file=paste(Project, "_rlogTransformed_proteinCodingOnly.tsv", sep=""), sep="\t")



rld.tissuecell <- rlogTransformation(dds.tissuecell, blind=T)
rld.tissue     <- rlogTransformation(dds.tissue, blind=T)


vsd.group      <- varianceStabilizingTransformation(dds.group,      blind=T)
vsd.tissuecell <- varianceStabilizingTransformation(dds.tissuecell, blind=T)


message("+-------------------------------------------------------------------------------")
message("+ Create some pHeatMaps")
message("+-------------------------------------------------------------------------------")


res_uterus_vs_liver.ann

# Most variable
topVarGenes            <- head(order(rowVars(assay(rld.tissue)),decreasing=TRUE),200)
mat                    <- assay(rld.group)[ topVarGenes, ]

# Favourites
mat <- assay(rld.group)[ favoritegenesA, ]

#mat <- res_uterus_vs_liver[ favoritegenesA, ]

mat                    <- mat - rowMeans(mat)
df                     <- as.data.frame(colData(rld.group)[,c("cell", "tissue")])
mat.df                 <- as.data.frame(mat)
mat.df$ensembl_gene_id <- rownames(mat.df)
mat.ann                <- merge(mat.df, ensEMBL2id, by="ensembl_gene_id")
rownames(mat.ann)      <- mat.ann$external_gene_name
#mat.ann                <- mat.ann[order(sapply(mat.ann$ensembl_gene_id, function(x) which(x == favoritegenesA))), ]
mat.ann                <- within(mat.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))
matrix                 <- as.matrix(mat.ann)


#ScaleCols <- colorRampPalette(colors = c("red3","white","royalblue4"))(255)
ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(255)
AnnotCols <- list( tissue=c("liver"="lightgrey", "uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"))

#pdf(paste(Project, "_DESeq_pHeatMap_rldgroup.pdf", sep=""), width=10,height=7, onefile=FALSE)
#par(bg=NA)
pheatmap(matrix, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, annotation_colors=AnnotCols, col=ScaleCols, main=paste(Project, " DESeq2 pHeatMap rld.tissue", sep=""))
#dev.off()

#
# By Log Fold Change
#
logFoldChanceCutOff <- 2.5

res_uterus_vs_liver.ann.bylfc     <- res_uterus_vs_liver.ann[ order( - abs(res_uterus_vs_liver.ann$log2FoldChange) ), ] 
res_uterus_vs_liver.ann.toplfc    <-  res_uterus_vs_liver.ann.bylfc[ !(is.na(res_uterus_vs_liver.ann.bylfc$log2FoldChange)) & 
                                                                   (abs(res_uterus_vs_liver.ann.bylfc$log2FoldChange) >= logFoldChanceCutOff) & 
                                                                  !(is.na(res_uterus_vs_liver.ann.bylfc$padj)) & 
                                                                   (res_uterus_vs_liver.ann.bylfc$padj <= significance), ] 

rows                    <- match(res_uterus_vs_liver.ann.toplfc$ensembl_gene_id, row.names(rld.tissue))
mat2                    <- assay(rld.tissue)[rows,]
df2                     <- as.data.frame(colData(rld.tissue)[,c("cell", "tissue")])
mat2.df                 <- as.data.frame(mat2)
mat2.df$ensembl_gene_id <- rownames(mat2.df)
mat2.ann                <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
rownames(mat2.ann)      <- mat2.ann$external_gene_name
mat2.ann                <- within(mat2.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))
matrix2                 <- as.matrix(mat2.ann)

ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(255)
AnnotCols <- list( tissue=c("liver"="lightgrey", "uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"))

pdf(paste0(Project, "_DESeq_pHeatMap_rld.tissue_lfc.", logFoldChanceCutOff, "_sig.", significance,  ".pdf"), width=10,height=7, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df2, col=ScaleCols, annotation_colors=AnnotCols, main=paste(Project, " DESeq2 pHeatMap rld.tissue log2FC=", logFoldChanceCutOff, " padj=", significance,  sep=""))
dev.off()




favourites <- c("ENSMUSG00000062960", "ENSMUSG00000029231", "ENSMUSG00000047414", "ENSMUSG00000006369", "ENSMUSG00000017737",
                "ENSMUSG00000035493", "ENSMUSG00000022150", "ENSMUSG00000039239", "ENSMUSG00000041324", "ENSMUSG00000052187",
                "ENSMUSG00000055609", "ENSMUSG00000052217", "ENSMUSG00000072115", "ENSMUSG00000004328", "ENSMUSG00000027993", 
                "ENSMUSG00000001507" )
#res_uterus_vs_liver.ann        <- functionGetSigDEG(as.data.frame(res_uterus_vs_liver), Project,"res_uterus_vs_liver",
#                                                    significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
#res_uterus_vs_liver.ann.bylfc  <- res_uterus_vs_liver.ann[ order( - abs(res_uterus_vs_liver.ann$log2FoldChange) ), ] 
#rownames(res_uterus_vs_liver.ann.bylfc) <- res_uterus_vs_liver.ann.bylfc$ensembl_gene_id
#res_uterus_vs_liver.ann.bylfc           <- res_uterus_vs_liver.ann.bylfc[favourites,]


res_uterus_vs_liver.df                 <- as.data.frame(res_uterus_vs_liver)
res_uterus_vs_liver.df                 <- res_uterus_vs_liver.df[ favourites, ]
res_uterus_vs_liver.df$ensembl_gene_id <- rownames(res_uterus_vs_liver.df)

rows                    <- match(res_uterus_vs_liver.df$ensembl_gene_id, row.names(rld.tissue))
mat2                    <- assay(rld.tissue)[rows,]
df2                     <- as.data.frame(colData(rld.tissue)[,c("cell", "tissue")])
df2.means               <- as.data.frame( list( "tissue"=c("uterus","liver")) )
rownames(df2.means)     <- c("uterus","liver")

mat2.df                 <- as.data.frame(mat2)
mat2.df$ensembl_gene_id <- rownames(mat2.df)
mat2.ann                <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
mat2.ann                <- mat2.ann[order(sapply(mat2.ann$ensembl_gene_id, function(x) which(x == favourites))), ]
rownames(mat2.ann)      <- mat2.ann$external_gene_name
mat2.ann                <- within(mat2.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))

mat2.ann.means         <- mat2.ann
mat2.ann.means$uterus  <- rowMeans(subset(mat2.ann, select = c("01","02","06","07","08","12","13")), na.rm = TRUE)
mat2.ann.means$liver   <- rowMeans(subset(mat2.ann, select = c("04","05","09","10","14","15")), na.rm = TRUE)
mat2.ann.means         <- mat2.ann.means[ , c("uterus","liver") ]

matrix2                 <- as.matrix(mat2.ann)
matrix2.means           <- as.matrix(mat2.ann.means)


ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(255)

ScaleCols <- colorRampPalette(colors = c("red","white","blue"))(255)

AnnotCols <- list( tissue=c("liver"="lightgrey", "uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"))

pdf(paste0(Project, "_DESeq_pHeatMap_rld.tissue_res_uterus_vs_liver_CustomList.pdf"), width=7.5,height=12.5, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 40, main=paste0(Project, " custom", "\nrld.tissue/res_uterus_vs_liver"))
dev.off()

pdf(paste0(Project, "_DESeq_pHeatMap_rld.tissue_res_uterus_vs_liver_CustomList_means.pdf"), width=5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2.means, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=FALSE, annotation_col=df2.means, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 40, main=paste0(Project, " custom", "\nrld.tissue/res_uterus_vs_liver"))
dev.off()



#
# Pairwise group comaprisons
#

logFoldChanceCutOff <- 5
significance        <- 0.01

res_uterus_cNK_vs_uterus_trNK.ann              <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_uterus_trNK),             Project,"res_uterus_cNK_vs_uterus_trNK",             significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_cNK_vs_uterus_ILC1.ann              <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_uterus_ILC1),             Project,"res_uterus_cNK_vs_uterus_ILC1",             significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_ILC1_vs_uterus_trNK.ann             <- functionGetSigDEG(as.data.frame(res_uterus_ILC1_vs_uterus_trNK),            Project,"res_uterus_ILC1_vs_uterus_trNK",            significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_liver_cNK_vs_liver_ILC1.ann                <- functionGetSigDEG(as.data.frame(res_liver_cNK_vs_liver_ILC1),               Project,"res_liver_cNK_vs_liver_ILC1",               significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_cNK_vs_liver_cNK.ann                <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_liver_cNK),               Project,"res_uterus_cNK_vs_liver_cNK",               significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_ILC1_vs_liver_ILC1.ann              <- functionGetSigDEG(as.data.frame(res_uterus_ILC1_vs_liver_ILC1),             Project,"res_uterus_ILC1_vs_liver_ILC1",             significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_trNK_vs_liver_ILC1.ann              <- functionGetSigDEG(as.data.frame(res_uterus_trNK_vs_liver_ILC1),             Project,"res_uterus_trNK_vs_liver_ILC1",             significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_trNK_vs_liver_cNK.ann               <- functionGetSigDEG(as.data.frame(res_uterus_trNK_vs_liver_cNK),              Project,"res_uterus_trNK_vs_liver_cNK",              significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_liver_cNK_vs_liver_ILC1.ann                <- functionGetSigDEG(as.data.frame(res_liver_cNK_vs_liver_ILC1),               Project,"res_liver_cNK_vs_liver_ILC1",               significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))
res_uterus_cNK_vs_liver_ILC1.ann               <- functionGetSigDEG(as.data.frame(res_uterus_cNK_vs_liver_ILC1),              Project,"res_uterus_cNK_vs_liver_ILC1",              significance,logFoldChanceCutOff,ensEMBL2id, as.data.frame(rld.group.ann))

genes2plot <- unique( c(res_uterus_cNK_vs_uterus_trNK.ann$ensembl_gene_id, 
                        res_uterus_cNK_vs_uterus_ILC1.ann$ensembl_gene_id,
                        res_uterus_ILC1_vs_uterus_trNK.ann$ensembl_gene_id,
                        res_liver_cNK_vs_liver_ILC1.ann$ensembl_gene_id,
                        res_uterus_cNK_vs_liver_cNK.ann$ensembl_gene_id,
                        res_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                        res_uterus_trNK_vs_liver_ILC1.ann$ensembl_gene_id,
                        res_uterus_trNK_vs_liver_cNK.ann$ensembl_gene_id,
                        res_liver_cNK_vs_liver_ILC1.ann$ensembl_gene_id,
                        res_uterus_cNK_vs_liver_ILC1.ann$ensembl_gene_id
) )
length(genes2plot)

rows3                   <- match(genes2plot, row.names(rld.group))
mat3                    <- assay(rld.group)[rows3,]
df3                     <- as.data.frame(colData(rld.group)[,c("cell", "tissue")])
mat3.df                 <- as.data.frame(mat3)
mat3.df$ensembl_gene_id <- rownames(mat3.df)
mat3.ann                <- merge(mat3.df, ensEMBL2id, by="ensembl_gene_id")
mat3.ann                <- mat3.ann[ !duplicated(mat3.ann[,c("external_gene_name")]) , ]
rownames(mat3.ann)      <- mat3.ann$external_gene_name
mat3.ann                <- within(mat3.ann, rm("ensembl_gene_id","external_gene_name","description","entrezgene","chromosome_name","gene_biotype","Length"))
matrix3                 <- as.matrix(mat3.ann)

colors<-colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(255)

pdf(paste0(Project, "_DESeq_pHeatMap_rld.group_pairwiseDEGs_lfc.", logFoldChanceCutOff, "_sig.", significance,  ".pdf"), width=10,height=7, onefile=FALSE)
par(bg=NA)
pheatmap(matrix3, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df3, 
         annotation_colors=list(tissue=c(liver="lightgrey",uterus="darkgrey"), cell=c(cNK="orange", ILC1="firebrick3", trNK="steelblue3")),
         cutree_cols=6, border_color="white",
         col=colors, main=paste(Project, " DESeq2 pHeatMap rld.groups log2FC=", logFoldChanceCutOff, " padj=", significance,  sep=""))
dev.off()








message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

rv     = rowVars(assay(rld.group))
select = order(rv, decreasing = TRUE)[seq_len(min(47000, length(rv)))]
pca    = prcomp(t(assay(rld.group)[select, ]))
fac    = factor(apply(as.data.frame(colData(rld.group)[, "group", drop = FALSE]), 1, paste, collapse = " : "))
                        
pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(sampleFiles=sampleTable$fileName, pca$x, sampleCell=sampleTable$cell, sampleName=sampleTable$sampleName, sampleTissue=sampleTable$tissue)

pdf(paste(Project, "_DESeq2_Annotated_PCA_All_nolabels.pdf", sep=""),width=10,height=7)
par(bg=NA)
  ggplot(data=scores, 
         aes(x = scores$PC1, y = scores$PC2, col = factor(scores$sampleCell), shape=(factor(scores$sampleTissue)), 
         label=scores$sampleFiles )) +
  geom_point(size = 5 ) + 
  xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste(Project, " PCA (All Genes)", sep="")) +
  scale_shape_manual(name="Tissue", values = c(17, 16)) + 
  scale_fill_manual(name="Cell Type",  values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +
  scale_colour_manual(name="Cell Type", values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +
  theme(text = element_text(size=elementTextSize)) 
dev.off()

pdf(paste(Project, "_DESeq2_Annotated_PCA_All_labels.pdf", sep=""),width=10,height=7)
par(bg=NA)
ggplot(data=scores, 
       aes(x = scores$PC1, y = scores$PC2, col = factor(scores$sampleCell), shape=(factor(scores$sampleTissue)), label=scores$sampleFiles )) +
  geom_point(size = 5 ) + #geom_text(aes(label=scores$sampleName), nudge_x=runif(1, 2, 4), nudge_y = runif(1, -2.5, 2.5), check_overlap = FALSE, size=3) +
  geom_text_repel(aes(label=scores$sampleName), box.padding = unit(1.0, "lines")) +
  xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste(Project, " PCA (All Genes)", sep="")) +
  scale_shape_discrete(name="Tissue") + scale_colour_manual(name="Cell Type", values = c("#1B9E77","#F27314", "purple", "red", "blue")) +
  theme(text = element_text(size=elementTextSize)) 
dev.off()

message("+-------------------------------------------------------------------------------")
message("+ Explore PCA Loadings")
message("+-------------------------------------------------------------------------------")


loadings                 <- as.data.frame(pca$rotation)
loadings$ensembl_gene_id <- rownames(loadings)
loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")



pca.1         <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <-  pca.1[c(1:25),]
pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.2         <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <-  pca.2[c(1:25),]
pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.3         <-  loadings[ order(loadings$PC3,decreasing=TRUE), ]
pca.3.25      <-  pca.3[c(1:25),]
pca.3.25.plot <- ggplot(data=pca.3.25, aes(x=factor(pca.3.25$external_gene_name,levels=unique(pca.3.25$external_gene_name)), y=PC3)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.4         <-  loadings[ order(loadings$PC4,decreasing=TRUE), ]
pca.4.25      <-  pca.4[c(1:25),]
pca.4.25.plot <- ggplot(data=pca.4.25, aes(x=factor(pca.4.25$external_gene_name,levels=unique(pca.4.25$external_gene_name)), y=PC4)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


pdf(paste(Project, "_DESeq2_PCA_Loadings.pdf", sep=""),width=10,height=7)
par(bg=NA)
plot_grid(pca.1.25.plot, pca.2.25.plot, pca.3.25.plot, pca.4.25.plot, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)
dev.off()

pc1var  <- round(summary(pca)$importance[2,1]*100,  digits=1)
pc2var  <- round(summary(pca)$importance[2,2]*100,  digits=1)
pc3var  <- round(summary(pca)$importance[2,3]*100,  digits=1)
pc4var  <- round(summary(pca)$importance[2,4]*100,  digits=1)
pc5var  <- round(summary(pca)$importance[2,5]*100,  digits=1)
pc6var  <- round(summary(pca)$importance[2,6]*100,  digits=1)
pc7var  <- round(summary(pca)$importance[2,7]*100,  digits=1)
pc8var  <- round(summary(pca)$importance[2,8]*100,  digits=1)
pc9var  <- round(summary(pca)$importance[2,9]*100,  digits=1)
pc10var <- round(summary(pca)$importance[2,10]*100, digits=1)

pcaVar = data.frame(pcs =c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),  
                    vals=c(pc1var, pc2var, pc3var, pc4var, pc5var, pc6var, pc7var, pc8var, pc9var, pc10var))


pdf(paste(Project, "_DESeq2_PCA_VarianceExplained.pdf", sep=""),width=10,height=7)
par(bg=NA)
ggplot(pcaVar, aes(factor(pcaVar$pcs,levels=unique(pcaVar$pcs)), vals)) + geom_col() + xlab("Principal Components") + ylab("Variance Explained") + ggtitle(paste0(Project, " Principal Components"))
dev.off()

library("scatterplot3d")


head(scores, 15)

scores$colours <- c("#1B9E77","#F27314","purple", "#1B9E77","#1B9E77","#F27314", "purple","purple","#1B9E77", "#F27314","purple","purple", "#1B9E77")

pdf(paste(Project, "_DESeq2_PCA_3D.pdf", sep=""),width=10,height=7)
par(bg=NA)
scatterplot3d(scores$PC1, scores$PC3, scores$PC2, color=scores$colours, highlight.3d = FALSE, col.axis = "black", col.grid = "grey", main = paste0(Project, " 3D PCA"), pch = 20, xlab = "PC1", ylab="PC3", zlab="PC2")
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ Create Per Gene Plots")
message("+-------------------------------------------------------------------------------")


makeGeneCountPlot <- function(dds,Project,gene2plot,outdir) {
  #
  # Plot the normalised read counts for a specified gene
  #
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(dds, gene=gene2plot, intgroup=c("group"), normalized=TRUE, returnData=TRUE)
  
  pdf(paste(outdir, Project, "-DGE_", gene2plot, "_collated.pdf", sep=""),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({
  ggplot(t2, aes(x=group, y=count, fill=group)) + geom_violin(trim=TRUE, alpha=.75)  + geom_boxplot(width = 0.1, fill='white') + 
  geom_point(position=position_jitter(w=0.1,h=0), alpha=0.25) +
  ggtitle(paste(Project, " ::: ", gene2plot, " ::: ", genename2plot, sep="")) + ylab(bquote('Normalised Read Count (' ~log[10]~ ')')) +
  #scale_x_discrete(limits=c("Virgin", "E9.5","E18.5")) + 
  scale_y_log10() +
  #scale_fill_manual(values = c("#1B9E77","#F27314", "purple")) + 
  theme(text = element_text(size=elementTextSize), legend.position="none") })
  dev.off()
  
  t2$samples    <- rownames(t2)
  
  pdf(paste(outdir, Project, "-DGE_", gene2plot, "_individual.pdf", sep=""),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ ggplot(t2, aes(x=samples, y=count, fill=group)) + geom_bar(stat="identity", alpha=.75) +
  #scale_fill_manual(values = c("#1B9E77","#F27314", "purple")) +
  ylab(bquote('Normalised Read Count')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=elementTextSize)) + ggtitle(paste(Project, " ::: ", gene2plot, " ::: ", genename2plot, sep="")) })
  dev.off()
}


makeGeneCountPlot(dds.group,Project,'ENSMUSG00000025044')
  
makeGeneCountPlot(dds.group,Project,'ENSMUSG00000028150')


message("+-------------------------------------------------------------------------------")
message("+ Create custom faceted selected gene plots (violin/dotplot")
message("+-------------------------------------------------------------------------------")


makeGeneSetExpressionPlot <- function(dds,Project,SelectedGeneList,SetName) {
  elementTextSize <- 8
  norm_counts                   <- dds
  SelectedGenes                 <- as.data.frame( ( norm_counts[ rownames(norm_counts) %in% SelectedGeneList,] ) )
  SelectedGenes                 <- SelectedGenes[SelectedGeneList, ]
  SelectedGenes$ensembl_gene_id <- rownames(SelectedGenes)
  SelectedGenes                 <- merge(SelectedGenes, ensEMBL2id, by="ensembl_gene_id")
  SelectedGenes$ensembl_gene_id <- ordered(SelectedGenes$ensembl_gene_id, SelectedGeneList ) 
  SelectedGenes                 <- SelectedGenes[with(SelectedGenes, order(ensembl_gene_id)),]
  rownames(SelectedGenes)       <- SelectedGenes$external_gene_name
  facet_order                   <-  SelectedGenes$external_gene_name
  SelectedGenes                 <- as.data.frame(t( within(SelectedGenes, rm("ensembl_gene_id","external_gene_name","description","entrezgene","chromosome_name","gene_biotype","Length"))))
  SelectedGenes$cell            <- sampleTable$cell
  SelectedGenes$tissue          <- sampleTable$tissue
  
  #SelectedGenes[SelectedGenes ==  0] <- NA
  SelectedGenes.mlt                  <- melt(SelectedGenes)
  SelectedGenes.mlt$variable         <- ordered(SelectedGenes.mlt$variable, facet_order ) 
  SelectedGenes.mlt                  <- SelectedGenes.mlt[with(SelectedGenes.mlt, order(variable)),]
  
  print({ head(SelectedGenes, 15)})
  
  plot <- ggplot(SelectedGenes.mlt, aes(x=cell, y=(log2(value+1)), colour=(log2(value+1)))) + geom_jitter(width = 0.15, size=2.5) + expand_limits( y = 0) +
    scale_colour_gradientn("log2(FPKM+1)", colours=c("yellow","red","black")) + facet_grid(variable ~ tissue, scales="free_x") +
    theme(axis.text.x = element_text(angle = 90, size=elementTextSize), axis.text.y = element_text(size=elementTextSize),
          strip.text.x = element_text(size=8, margin = margin(.1, 0, .1, 0, "cm")), strip.text.y = element_text(size=8, margin = margin(0, .1, 0, .1, "cm")),
          panel.grid.major = element_line(colour = "grey10", size=0.1)) +
    xlab("") + ylab("log2(FPKM+1)") + ylim(0,10)

  myHeight <- (1.3*nrow(as.data.frame(SelectedGeneList)))
  message(myHeight)
  
  pdf(paste(Project, "-DGE_GeneSetExpressionPlot_", SetName, ".pdf", sep=""),width=3.5,height=myHeight, onefile=FALSE)
  par(bg=NA)
  print({ plot })
  dev.off()
  
  return(plot)
}


selected.TF    <- c("ENSMUSG00000028150", "ENSMUSG00000019256", "ENSMUSG00000001444", "ENSMUSG00000032446", "ENSMUSG00000015619", "ENSMUSG00000032238")
plt.tf <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.TF,"TF")
plt.tf



selected.SRT   <- c("ENSMUSG00000042284", "ENSMUSG00000030724", "ENSMUSG00000032093", "ENSMUSG00000030325", "ENSMUSG00000062524",
                    "ENSMUSG00000030361", "ENSMUSG00000023274", "ENSMUSG00000024669", "ENSMUSG00000053977")
plt.SRT <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.SRT,"Sorting")

selected.TR   <- c("ENSMUSG00000037940", "ENSMUSG00000048521", "ENSMUSG00000019960", "ENSMUSG00000026573", "ENSMUSG00000074063",
                   "ENSMUSG00000018381", "ENSMUSG00000049410", "ENSMUSG00000021591", "ENSMUSG00000045087", "ENSMUSG00000040596",
                   "ENSMUSG00000057191", "ENSMUSG00000045382", "ENSMUSG00000044199", "ENSMUSG00000035900", "ENSMUSG00000042978",
                   "ENSMUSG00000055148", "ENSMUSG00000110753", "ENSMUSG00000032412", "ENSMUSG00000022696", "ENSMUSG00000000782",
                   "ENSMUSG00000043008", "ENSMUSG00000085603", "ENSMUSG00000045092", "ENSMUSG00000027330", "ENSMUSG00000036006",
                   "ENSMUSG00000031453", "ENSMUSG00000032946",  "ENSMUSG00000020901", "ENSMUSG00000004568")
plt.tr <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.TR,"TR")

selected.MY   <- c("ENSMUSG00000022582", "ENSMUSG00000019966", "ENSMUSG00000005339", "ENSMUSG00000039013", "ENSMUSG00000030786",
                   "ENSMUSG00000030789", "ENSMUSG00000019987", "ENSMUSG00000014361", "ENSMUSG00000004730", "ENSMUSG00000026395",
                   "ENSMUSG00000051457", "ENSMUSG00000034783", "ENSMUSG00000079018", "ENSMUSG00000075122", "ENSMUSG00000022901",
                   "ENSMUSG00000026712", "ENSMUSG00000031494", "ENSMUSG00000065987", "ENSMUSG00000040165", "ENSMUSG00000031495",
                   "ENSMUSG00000040197", "ENSMUSG00000051906", "ENSMUSG00000014773", "ENSMUSG00000026923", "ENSMUSG00000035385",
                   "ENSMUSG00000008845", "ENSMUSG00000024621", "ENSMUSG00000059326", "ENSMUSG00000071713", "ENSMUSG00000028859")
plt.MY <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.MY,"Myeloid")

selected.TLR <- c("ENSMUSG00000044827", "ENSMUSG00000027995", "ENSMUSG00000031639", "ENSMUSG00000039005", "ENSMUSG00000079164",
                  "ENSMUSG00000051498", "ENSMUSG00000044583", "ENSMUSG00000040522", "ENSMUSG00000045322", "ENSMUSG00000051969",
                  "ENSMUSG00000062545", "ENSMUSG00000033777")
plt.TLR <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.TLR,"TLR")


selected.MHC <- c("ENSMUSG00000036594", "ENSMUSG00000073421", "ENSMUSG00000079547", "ENSMUSG00000037649", "ENSMUSG00000060586",
                  "ENSMUSG00000067341")
plt.MHC <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.MHC,"MHC")


selected.ILC <- c("ENSMUSG00000026069", "ENSMUSG00000034117", "ENSMUSG00000032011", "ENSMUSG00000003882")
plt.ILC <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.ILC,"ILC")

#selected.INT <- c("ENSMUSG00000055170", "ENSMUSG00000031712", "ENSMUSG00000020383", "ENSMUSG00000000869", "ENSMUSG00000036117",
#                  "ENSMUSG00000025746", "ENSMUSG00000074695", "ENSMUSG00000025383", "ENSMUSG00000025929", "ENSMUSG00000024578",
#                  "ENSMUSG00000046108", "ENSMUSG00000050222", "ENSMUSG00000025929", "ENSMUSG00000041872", "ENSMUSG00000016529",
#                  "ENSMUSG00000002603", "ENSMUSG00000027776", "ENSMUSG00000004296", "ENSMUSG00000018916", "ENSMUSG00000014599",
#                  "ENSMUSG00000038067", "ENSMUSG00000031750")
#plt.INT <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.INT,"INT")


message("+-------------------------------------------------------------------------------")
message("+ Create custom faceted selected gene plots (violin/dotplot")
message("+-------------------------------------------------------------------------------")


all.points.df       <- read.table("DataTables/180410_all_points_uterus.txt", header=TRUE, stringsAsFactors=FALSE)
all.points.df$label <- paste(all.points.df$cell, all.points.df$time, sep="|")
all.points.df       <- all.points.df[ mixedorder(all.points.df$label),]
print(all.points.df)

all.points.df.summary        <- ddply(all.points.df, c( "label"), summarise, mean = median(value))
all.points.df.summary$cell   <- c("cNK","cNK","cNK","cNK","cNK", "cNK","cNK","cNK","cNK", "cNK",
                                  "ILC1", "ILC1","ILC1","ILC1","ILC1","ILC1", "ILC1","ILC1","ILC1", "ILC1",
                                 "trNK", "trNK","trNK","trNK","trNK","trNK","trNK","trNK","trNK", "trNK")
all.points.df.summary        <- all.points.df.summary[ mixedorder(all.points.df.summary$label),]
all.points.df.summary$time   <- gsub(".*\\|", "", all.points.df.summary$label)
all.points.df.summary$time   <- factor(all.points.df.summary$time, levels=c( "3w", "5w", "8w", "5.5", "9.5", "13.5", "17.5", "1d","9d", "18d"))
all.points.df.summary        <- all.points.df.summary[order(all.points.df.summary$time),]
print(all.points.df.summary)

all.points.plot <- ggplot() +  
                   geom_path(data=all.points.df.summary, aes(x=factor(time), y=mean, group=cell, colour=cell))  +
                   geom_point(data=all.points.df.summary, aes(x=factor(time), y=mean, group=cell, colour=cell, size=mean)) + 
                   scale_size_area(max_size = 20) +
                   xlab("Time Point") + ylim(0,100) + ylab("Percentage") +
                   scale_x_discrete( limits = c("3w", "5w", "8w","5.5", "9.5", "13.5", "17.5", "1d","9d", "18d"),
                                     labels = c("3w", "5w", "8w","gd5.5", "gd9.5", "gd13.5", "gd17.5","1d-pp", "10d-pp", "18d-pp")) +
                   scale_colour_manual("Cell Type", values=c("#FFCC00", "#EA5160", "#0099FF")) +
                   theme_minimal(base_size=12) + guides(size=FALSE)

all.points.plot.indi <- all.points.plot + geom_point(data=all.points.df, aes(y = value, x = time, size=1, colour=cell), alpha=1.0, fill='white', shape=23) 


CairoPDF(file = paste0(Project, "-BubbleTimePlot-All.pdf"), width = 12, height = 6, onefile = TRUE, family = "Arial" )
par(bg=NA)
all.points.plot
dev.off()

CairoPDF(file = paste0(Project, "-BubbleTimePlot-All-withindividualpoints.pdf"), width = 12, height = 6, onefile = TRUE, family = "Arial" )
par(bg=NA)
all.points.plot.indi
dev.off()



inset.df       <- read.table("DataTables/180410_nBF.txt", header=TRUE, stringsAsFactors=FALSE)
inset.df$label <- paste(inset.df$cell, inset.df$time, sep="|")
inset.df       <- inset.df[ mixedorder(inset.df$label),]
print(inset.df)

inset.df.summary        <- ddply(inset.df, c( "label"), summarise, mean = median(value))
inset.df.summary$cell   <- c("cNK", "ILC1", "trNK")
inset.df.summary        <- inset.df.summary[ mixedorder(inset.df.summary$label),]
inset.df.summary$time   <- gsub(".*\\|", "", inset.df.summary$label)
inset.df.summary$time   <- factor(inset.df.summary$time, levels=c( "10d"))
inset.df.summary        <- inset.df.summary[order(inset.df.summary$time),]
print(inset.df.summary)

inset.df.plot <- ggplot() +  
                 geom_path(data=inset.df.summary, aes(x=factor(time), y=mean, group=cell, colour=cell))  + 
                 geom_point(data=inset.df.summary, aes(x=factor(time), y=mean, group=cell, colour=cell, size=mean)) + scale_size_area(max_size = 20) +
                 xlab("Time Point") + ylim(0,100) + ylab("Percentage") + 
                 scale_x_discrete( limits = c("10d"),
                                   labels = c("10d-pp(nBF)")) +
                 scale_colour_manual("Cell Type", values=c("#FFCC00", "#EA5160", "#0099FF")) +
                 theme_minimal(base_size=22) + guides(size=FALSE)

inset.df.plot.indi <- inset.df.plot + geom_point(data=inset.df, aes(y = value, x = time, size=1, colour=cell), alpha=1.0, fill='white', shape=23) 

CairoPDF(file = paste0(Project, "-BubbleTimePlot-Inset.pdf"), width = 4, height = 6, onefile = TRUE, family = "Arial" )
par(bg=NA)
inset.df.plot
dev.off()

CairoPDF(file = paste0(Project, "-BubbleTimePlot-Inset-withindividualpoints.pdf"), width = 4, height = 6, onefile = TRUE, family = "Arial" )
par(bg=NA)
inset.df.plot.indi
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ CIBERSORT + ImmuCC")
message("+-------------------------------------------------------------------------------")

raw.counts                    <- as.data.frame( counts(ddsHTSeq.tissuecell, normalized = FALSE, replaced = FALSE) )
raw.counts$ensembl_gene_id    <- rownames(raw.counts)
raw.counts.ann                <- merge(raw.counts,ensEMBL2id, by="ensembl_gene_id" )
raw.counts.ann                <- raw.counts.ann[,c(15,2:14)]
colnames(raw.counts.ann)      <- c("Gene",  "cNKuterus.01", "trNKuterus.02", "ILC1liver.04", "cNKliver.05", "cNKuterus.06", 
                                   "trNKuterus.07", "ILC1uterus.08", "ILC1liver.09", "cNKliver.10", "trNKuterus.12", 
                                   "ILC1uterus.13","ILC1liver.14", "cNKliver.15")
raw.counts.ann.uniq           <- raw.counts.ann[!duplicated(raw.counts.ann$Gene),]
rownames(raw.counts.ann.uniq) <- raw.counts.ann.uniq$Gene
raw.counts.ann.uniq           <- raw.counts.ann.uniq[,-c(1)]
head(raw.counts.ann.uniq)

write.csv(raw.counts.ann.uniq, file=paste(Project, "RawReadCounts_tissuecell.csv", sep=""))


message("+-------------------------------------------------------------------------------")
message("+ Deconvolve the RNA counts to immune cell type proportions")
message("+-------------------------------------------------------------------------------")

library(DeconRNASeq)

# Prep the ImmuCC data matrix for DeconRNASeq
Immune_cells_expr_matrix                <- read.csv("/Users/rhamilto/Documents/CTR-Data/ImmuCC/srep40508-s1.csv")
names(Immune_cells_expr_matrix)[1]      <- paste("external_gene_name") 
Immune_cells_expr_matrix_anno           <- unique(merge(Immune_cells_expr_matrix, ensEMBL2id, by = "external_gene_name"))
Immune_cells_expr_matrix_anno           <- Immune_cells_expr_matrix_anno[order(Immune_cells_expr_matrix_anno$entrezgene),]
Immune_cells_expr_matrix_anno           <- Immune_cells_expr_matrix_anno[!duplicated(Immune_cells_expr_matrix_anno$external_gene_name),]
rownames(Immune_cells_expr_matrix_anno) <- Immune_cells_expr_matrix_anno$ensembl_gene_id
Immune_cells_expr_matrix_anno           <- Immune_cells_expr_matrix_anno[,-c(1,27:31)]
head(Immune_cells_expr_matrix_anno)

# Prep the projects nomalised counts
NormCounts                        <- normCounts.tissuecell
NormCounts_Immune                 <- as.matrix(NormCounts[rownames(NormCounts) %in% rownames(Immune_cells_expr_matrix_anno),])
Immune_cells_expr_matrix_anno2    <- as.matrix(Immune_cells_expr_matrix_anno[rownames(Immune_cells_expr_matrix_anno) %in% rownames(NormCounts_Immune),])

head(NormCounts_Immune)
head(Immune_cells_expr_matrix_anno2)


makeCellTypeDeconvolutionPlot <- function(Cell_Type, Decon_results_df) {
  
  t2   <- subset(Decon_results_df, select=c(Cell_Type, "cell", "tissue", "group"))
  t2   <- t2[order(t2$cell),] 
  
  plot <-  ggplot(t2, aes(x=group, y=t2[[1]], fill=tissue)) + 
          # geom_violin(trim=TRUE, alpha=.5)  + 
           geom_boxplot(width = 0.1, alpha=0.5) + 
           geom_point(position=position_jitter(w=0.05,h=0), alpha=0.99, colour="red") +
           ggtitle(Cell_Type) +
           scale_fill_manual(values = c("purple3", "olivedrab2") )  + 
           labs(y = "Cell type fraction", x = "") +
           ylim(0,0.75) +
          # expand_limits(y = 0) +
           theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")

 return(plot)
}


# Run the DeconRNASeq algorithm
Decon_results    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(1:10)]), 
                                checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)

Decon_results_df1           <- as.data.frame(Decon_results$out.all)
rownames(Decon_results_df1) <- paste0(sampleTable$sampleName, "_", sampleTable$group)

Decon_results_df1$sampleName <- sampleTable$sampleName
Decon_results_df1$cell       <- sampleTable$cell
Decon_results_df1$tissue     <- sampleTable$tissue
Decon_results_df1$group      <- sampleTable$group






plot.Mast.Cells          <- makeCellTypeDeconvolutionPlot('Mast.Cells',          Decon_results_df1)
plot.Neutrophil.Cells    <- makeCellTypeDeconvolutionPlot('Neutrophil.Cells',    Decon_results_df1)
plot.Eosinophil.Cells    <- makeCellTypeDeconvolutionPlot('Eosinophil.Cells',    Decon_results_df1)
plot.B.Cells.Memory      <- makeCellTypeDeconvolutionPlot('B.Cells.Memory',      Decon_results_df1)
plot.B.Cells.Naive       <- makeCellTypeDeconvolutionPlot('B.Cells.Naive',       Decon_results_df1)
plot.Plasma.Cells        <- makeCellTypeDeconvolutionPlot('Plasma.Cells',        Decon_results_df1)
plot.T.Cells.CD8.Actived <- makeCellTypeDeconvolutionPlot('T.Cells.CD8.Actived', Decon_results_df1)
plot.T.Cells.CD8.Naive   <- makeCellTypeDeconvolutionPlot('T.Cells.CD8.Naive',   Decon_results_df1)
plot.T.Cells.CD8.Memory  <- makeCellTypeDeconvolutionPlot('T.Cells.CD8.Memory',  Decon_results_df1)
plot.M0.Macrophage       <- makeCellTypeDeconvolutionPlot('M0.Macrophage',       Decon_results_df1)


# Run the DeconRNASeq algorithm
Decon_results    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(11:20)]), 
                                checksig = T, known.prop = FALSE, use.scale = TRUE, fig = TRUE)

Decon_results_df2           <- as.data.frame(Decon_results$out.all)
rownames(Decon_results_df2) <- paste0(sampleTable$sampleName, "_", sampleTable$group)
Decon_results_df2$cell      <- sampleTable$cell
Decon_results_df2$tissue    <- sampleTable$tissue
Decon_results_df2$group     <- sampleTable$group


plot.M1.Macrophage          <- makeCellTypeDeconvolutionPlot('M1.Macrophage',          Decon_results_df2)
plot.M2.Macrophage          <- makeCellTypeDeconvolutionPlot('M2.Macrophage',          Decon_results_df2)
plot.Treg.Cells             <- makeCellTypeDeconvolutionPlot('Treg.Cells',             Decon_results_df2)
plot.T.Cells.CD4.Memory     <- makeCellTypeDeconvolutionPlot('T.Cells.CD4.Memory',     Decon_results_df2)
plot.T.Cells.CD4.Naive      <- makeCellTypeDeconvolutionPlot('T.Cells.CD4.Naive',      Decon_results_df2)
plot.T.Cells.CD4.Follicular <- makeCellTypeDeconvolutionPlot('T.Cells.CD4.Follicular', Decon_results_df2)
plot.Th1.Cells              <- makeCellTypeDeconvolutionPlot('Th1.Cells',              Decon_results_df2)
plot.Th17.Cells             <- makeCellTypeDeconvolutionPlot('Th17.Cells',             Decon_results_df2)
plot.Th2.Cells              <- makeCellTypeDeconvolutionPlot('Th2.Cells',              Decon_results_df2)
plot.Monocyte               <- makeCellTypeDeconvolutionPlot('Monocyte',               Decon_results_df2)

# Run the DeconRNASeq algorithm
Decon_results    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(21:25)]), 
                                checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)

Decon_results_df3           <- as.data.frame(Decon_results$out.all)
rownames(Decon_results_df3) <- paste0(sampleTable$sampleName, "_", sampleTable$group)
Decon_results_df3$cell      <- sampleTable$cell
Decon_results_df3$tissue    <- sampleTable$tissue
Decon_results_df3$group     <- sampleTable$group

plot.GammaDelta.T.Cells     <- makeCellTypeDeconvolutionPlot('GammaDelta.T.Cells', Decon_results_df3)
plot.NK.Resting             <- makeCellTypeDeconvolutionPlot('NK.Resting',         Decon_results_df3)
plot.NK.Actived             <- makeCellTypeDeconvolutionPlot('NK.Actived',         Decon_results_df3)
plot.DC.Actived             <- makeCellTypeDeconvolutionPlot('DC.Actived',         Decon_results_df3)
plot.DC.Immature            <- makeCellTypeDeconvolutionPlot('DC.Immature',        Decon_results_df3)


pdf(paste0(Project, "-ImmuCC_DeconRNASEq.pdf"),width=20,height=20, onefile=FALSE)
par(bg=NA)
plot_grid(plot.Mast.Cells, plot.Neutrophil.Cells, plot.Eosinophil.Cells, plot.B.Cells.Memory, plot.B.Cells.Naive,
          plot.Plasma.Cells, plot.T.Cells.CD8.Actived, plot.T.Cells.CD8.Naive, plot.T.Cells.CD8.Memory,plot.M0.Macrophage,
          plot.M1.Macrophage, plot.M2.Macrophage, plot.Treg.Cells, plot.T.Cells.CD4.Memory, plot.T.Cells.CD4.Naive,
          plot.T.Cells.CD4.Follicular, plot.Th1.Cells, plot.Th17.Cells, plot.Th2.Cells, plot.Monocyte,
          plot.GammaDelta.T.Cells, plot.NK.Resting, plot.NK.Actived, plot.DC.Actived, plot.DC.Immature,
          ncol = 5, nrow = 5)
dev.off()


head(Decon_results_df1)
head(Decon_results_df2)
head(Decon_results_df3)

Decon_results_df <- merge(Decon_results_df1, Decon_results_df2,by="row.names")
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL

Decon_results_df <- merge(Decon_results_df,  Decon_results_df3,by="row.names",all.x=TRUE)
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL


Decon_results_df.means      <- aggregate(Decon_results_df[, c(1:10,15:28)], list(Decon_results_df$group), mean)
Decon_results_df.means.melt <- melt(Decon_results_df.means)

Decon_results_df.vars      <- aggregate(Decon_results_df[, c(1:10,15:28)], list(Decon_results_df$group), var)
Decon_results_df.vars.melt <- melt(Decon_results_df.vars)

Decon_results_df.means.melt$vars      <- Decon_results_df.vars.melt$value

Decon_results_df.means.melt$normvars  <- Decon_results_df.means.melt$vars/max(Decon_results_df.means.melt$vars)

colnames(Decon_results_df.means.melt) <- c("Group", "Immune.Cell", "mean", "variance", "normalisedvariance")

write.csv(Decon_results_df.means.melt, file=paste0(Project, "_Decon_results_df.means.melt.csv"))




deconRNASeq.summary.table <- read.table("CTR_fc287_0001_Decon_results_df.means.melt.csv", sep=",", header=TRUE, row.names=1, stringsAsFactors = FALSE)

selectedOrder <- rev(c("Treg.Cells","Th2.Cells","Th17.Cells","Th1.Cells","T.Cells.CD8.Naive",
                     "T.Cells.CD8.Actived","T.Cells.CD8.Memory","T.Cells.CD4.Naive","T.Cells.CD4.Memory","T.Cells.CD4.Follicular",
                     "GammaDelta.T.Cells","NK.Resting","NK.Actived","M0.Macrophage","M1.Macrophage",
                     "M2.Macrophage","DC.Actived", "Monocyte","Mast.Cells","Neutrophil.Cells",
                     "Eosinophil.Cells","B.Cells.Naive","B.Cells.Memory","Plasma.Cells" ))

deconRNASeq.summary.table$Immune.Cell <- factor(deconRNASeq.summary.table$Immune.Cell, levels = selectedOrder )
deconRNASeq.summary.table$Immune.Cell <- gsub("\\.",  " ",      deconRNASeq.summary.table$Immune.Cell)
deconRNASeq.summary.table$Group       <- gsub("cNK",  "cNK\n",  deconRNASeq.summary.table$Group)
deconRNASeq.summary.table$Group       <- gsub("ILC1", "ILC1\n", deconRNASeq.summary.table$Group)
deconRNASeq.summary.table$Group       <- gsub("trNK", "trNK\n", deconRNASeq.summary.table$Group)

ggplot(deconRNASeq.summary.table, aes(y = Immune.Cell,x = Group)) +
  geom_point(aes(size = mean, colour = Group, alpha=(1-normalisedvariance))) + 
  xlab("") + ylab("") +
  guides(alpha=guide_legend(title="1-(Normalised Variance)")) +
  guides(size=guide_legend(title="Mean Proportion")) +
  guides(colour=guide_legend(title="Group")) +  
  guides(colour=FALSE) +
  theme_bw()



message("+-------------------------------------------------------------------------------")
message("+ Analysis Complete")
message("+-------------------------------------------------------------------------------")


#------------------------------------------------------------------------------
# FIN
#------------------------------------------------------------------------------

