#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# RNA-Seq Analysis to accompany:
#
# Filipovic et al, 2018      
# Molecular definition of group 1 innate lymphoid cells in the mouse uterus
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2018-Filipovic-Colucci
#
# CTR Code: CTR_fc287_0001
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics 
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------


#
# initial install of packages
#
#source("http://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
#biocLite("apeglm")
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
library("GenomicRanges")
library("rtracklayer")
library("ggdendro")
library("ggalt")
library("apeglm")


message("+-------------------------------------------------------------------------------")
message("+ Set up some constants e.g. base directories")
message("+-------------------------------------------------------------------------------")

Project         <- "2018-Filipovic-Colucci"
Base.dir        <- getwd() 
setwd(Base.dir)
HTSeq.dir       <- paste(Base.dir,"/HTSeq_Counts", sep="")

# Set GTFfile to point to the file on your local computer e.g. 
GTFfile         <- "/usr/local/Genomes/Mus_musculus/GRCm38/Mus_musculus.GRCm38.84.gtf"

elementTextSize <- 10

significance    <- 0.05
foldchange      <- 2


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


message("+-------------------------------------------------------------------------------")
message("+ Retrieve average transcript lengths")
message("+-------------------------------------------------------------------------------")
# From: http://seqanswers.com/forums/archive/index.php/t-39797.html

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





message("+-------------------------------------------------------------------------------")
message("+ Perform Comparisons of defined groups, tissues, cells")
message("+-------------------------------------------------------------------------------")

DESeq2Version <- "1.18.0"

res_uterus_vs_liver                                                                    <- results(dds.tissuecell,   contrast=c("tissue", "uterus",  "liver"))
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_vs_liver                        <- lfcShrink(dds.tissuecell, contrast=c("tissue", "uterus",  "liver"),        type="normal") }


res_uterus_cNK_vs_uterus_trNK                                                          <- results(dds.group,        contrast=c("group", "cNKuterus", "trNKuterus")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_cNK_vs_uterus_trNK              <- lfcShrink(dds.group,      contrast=c("group", "cNKuterus", "trNKuterus"),  type="normal") }

res_uterus_trNK_vs_uterus_cNK                                                          <- results(dds.group,        contrast=c("group", "trNKuterus", "cNKuterus")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_trNK_vs_uterus_cNK              <- lfcShrink(dds.group,      contrast=c("group", "trNKuterus", "cNKuterus"),  type="normal") }

res_uterus_cNK_vs_uterus_ILC1                                                          <- results(dds.group,        contrast=c("group", "cNKuterus", "ILC1uterus")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_cNK_vs_uterus_ILC1              <- lfcShrink(dds.group,      contrast=c("group", "cNKuterus", "ILC1uterus"),  type="normal") }

res_uterus_ILC1_vs_uterus_trNK                                                         <- results(dds.group,        contrast=c("group", "ILC1uterus", "trNKuterus")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_ILC1_vs_uterus_trNK             <- lfcShrink(dds.group,      contrast=c("group", "ILC1uterus", "trNKuterus"), type="normal") }

res_liver_cNK_vs_liver_ILC1                                                            <- results(dds.group,        contrast=c("group", "cNKliver", "ILC1liver")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_liver_cNK_vs_liver_ILC1                <- lfcShrink(dds.group,      contrast=c("group", "cNKliver", "ILC1liver"),    type="normal") }

res_uterus_cNK_vs_liver_cNK                                                            <- results(dds.group,        contrast=c("group", "cNKuterus", "cNKliver")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_cNK_vs_liver_cNK                <- lfcShrink(dds.group,      contrast=c("group", "cNKuterus", "cNKliver"),    type="normal") }

res_uterus_ILC1_vs_liver_ILC1                                                          <- results(dds.group,        contrast=c("group", "ILC1uterus", "ILC1liver")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_ILC1_vs_liver_ILC1              <- lfcShrink(dds.group,      contrast=c("group", "ILC1uterus", "ILC1liver"),  type="normal") }

res_uterus_trNK_vs_liver_ILC1                                                          <- results(dds.group,        contrast=c("group", "trNKuterus", "ILC1liver")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_trNK_vs_liver_ILC1              <- lfcShrink(dds.group,      contrast=c("group", "trNKuterus", "ILC1liver"),  type="normal") }

res_uterus_trNK_vs_liver_cNK                                                           <- results(dds.group,        contrast=c("group", "trNKuterus", "cNKliver")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_trNK_vs_liver_cNK               <- lfcShrink(dds.group,      contrast=c("group", "trNKuterus", "cNKliver"),   type="normal") }

res_liver_cNK_vs_liver_ILC1                                                            <- results(dds.group,        contrast=c("group", "cNKliver", "ILC1liver"))
if(packageVersion("DESeq2")>DESeq2Version){ res_liver_cNK_vs_liver_ILC1                <- lfcShrink(dds.group,      contrast=c("group", "cNKliver", "ILC1liver"),    type="normal") }

res_uterus_cNK_vs_liver_ILC1                                                           <- results(dds.group,        contrast=c("group", "cNKuterus", "ILC1liver")) 
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_cNK_vs_liver_ILC1               <- lfcShrink(dds.group,      contrast=c("group", "cNKuterus", "ILC1liver"),   type="normal") }


res_uterus_trNK_uterus_ILC1_vs_liver_ILC1                                              <- results(dds.collated1,    contrast=c("collated1", "utrNK_uILC1", "lILC1"))
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_trNK_uterus_ILC1_vs_liver_ILC1  <- lfcShrink(dds.collated1,  contrast=c("collated1", "utrNK_uILC1", "lILC1"), type="normal") }

res_uterus_ILC1_uterus_trNK_vs_uterus_cNK                                              <- results(dds.collated2,    contrast=c("collated2", "uILC1_utrNK", "ucNK"))
if(packageVersion("DESeq2")>DESeq2Version){ res_uterus_ILC1_uterus_trNK_vs_uterus_cNK  <- lfcShrink(dds.collated2,  contrast=c("collated2", "uILC1_utrNK", "ucNK"),  type="normal") }



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



message("+-------------------------------------------------------------------------------")
message("+ Create MA Plots")
message("+-------------------------------------------------------------------------------")

favoritegenesA <- NULL
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000031375") #Bgn
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000025473") #Adam8
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000035385") #Ccl2
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000056529") #Ptafr
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000036006") #Fam65b
favoritegenesA <- append(favoritegenesA, "ENSMUSG00000027009") #Itga4

favoritegenesB <- NULL
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000027765") #P2ry1   
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000029470") #P2rx4
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000029832") #Nfe2l3
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000015839") #Nfe2l2
favoritegenesB <- append(favoritegenesB, "ENSMUSG00000025993") #Slc40a1

favoritegenesC <- NULL
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000021423") #Ly86
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000044827") #Tlr1
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000027995") #Tlr2
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000040522") #Tlr8
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000045322") #Tlr9
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000024164") #C3
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000031443") #F7
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000031444") #F10
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000025314") #Ptprj
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000040552") #C3ar1
favoritegenesC <- append(favoritegenesC, "ENSMUSG00000049130") #C5ar1

favoritegenesD <- NULL
favoritegenesD <- append(favoritegenesD, "ENSMUSG00000021822") #Plau
favoritegenesD <- append(favoritegenesD, "ENSMUSG00000025856") #Pdgfa
favoritegenesD <- append(favoritegenesD, "ENSMUSG00000027200") #Sema6D
favoritegenesD <- append(favoritegenesD, "ENSMUSG00000020689") #Itgb3
favoritegenesD <- append(favoritegenesD, "ENSMUSG00000021796") #Bmpr1a

megafavorite <- NULL
megafavorite <- append(megafavorite, favoritegenesA)
megafavorite <- append(megafavorite, favoritegenesB)
megafavorite <- append(megafavorite, favoritegenesC)
megafavorite <- append(megafavorite, favoritegenesD)


functionCustomMAPlot.1colour <- function(results, Project, FigureID, Title, significance, log2FoldChange, xfavoritegenes) {
  # Get annotation infor for a favorite set of ensEMBL ids
  labeldata.rows            <- match(xfavoritegenes, row.names(results))
  labeldata                 <- results[labeldata.rows,]
  labeldata$ensembl_gene_id <- rownames(labeldata)
  labeldata.ann             <- merge(labeldata, ensEMBL2id, by="ensembl_gene_id")
  
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
    
    geom_label_repel(data=labeldata.ann,
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     fill='gray', colour='black', point.padding = unit(0.25, "lines"),  size=6, segment.size = 1, segment.color = 'darkred',  nudge_x = 0, nudge_y=0) +
    scale_x_log10() +
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") + ggtitle(paste(Project, " ", FigureID, "\n", Title, " [fc ", log2FoldChange, ", sig ", significance, "]", sep=""))
  
  pdf(paste(Project, "_", FigureID, ".pdf", sep=""),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt })
  dev.off()
  return (plt)
}

SuppFig.S3A <- functionCustomMAPlot.1colour(as.data.frame(res_uterus_trNK_vs_uterus_cNK), Project, "SuppFig.S3A", "res_uterus_trNK_vs_uterus_cNK", significance, foldchange, megafavorite)




message("+-------------------------------------------------------------------------------")
message("+ UpSetR")
message("+-------------------------------------------------------------------------------")


listInput <- list("uterus vs liver"                         = res_uterus_vs_liver.ann$ensembl_gene_id,
                  'uterine trNK and ILC1 vs liver ILC1'     = res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine trNK and ILC1 vs uterine cNK'    = res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.ann$ensembl_gene_id,
                  'uterine trNK vs uterine cNK'             = res_uterus_cNK_vs_uterus_trNK.ann$ensembl_gene_id,
                  'uterine ILC1 vs uterine trNK'            = res_uterus_ILC1_vs_uterus_trNK.ann$ensembl_gene_id,
                  'uterine ILC1 vs uterine cNK'             = res_uterus_cNK_vs_uterus_ILC1.ann$ensembl_gene_id,
                  'liver ILC1 vs liver cNK'                 = res_liver_cNK_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine cNK vs liver cNK'                = res_uterus_cNK_vs_liver_cNK.ann$ensembl_gene_id,
                  'uterine ILC1 vs liver ILC1'              = res_uterus_ILC1_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine trNK vs liver ILC1'              = res_uterus_trNK_vs_liver_ILC1.ann$ensembl_gene_id,
                  'uterine trNK vs liver cNK'               = res_uterus_trNK_vs_liver_cNK.ann$ensembl_gene_id
                )

upset(fromList(listInput), nsets = 11, sets = rev(colnames(fromList(listInput))), cutoff = 7,
      keep.order = TRUE, empty.intersections = "on", nintersects = 11, 
      sets.x.label = "Number of differentially expressed genes", mainbar.y.label = "Intersections of Differentially Expressed Genes")




pdf(paste(Project, "_Fig.3A_11sets.pdf", sep=""),width=10,height=7, onefile=FALSE)
par(bg=NA)
upset(fromList(listInput), nsets = 11, sets = rev(colnames(fromList(listInput))), cutoff = 7,
      keep.order = TRUE, empty.intersections = "on", nintersects = 11, 
      sets.x.label = "Number of differentially expressed genes", mainbar.y.label = "Intersections of Differentially Expressed Genes")
dev.off()


pdf(paste(Project, "_Fig.3A_5sets.pdf", sep=""),width=10,height=7, onefile=FALSE)
par(bg=NA)
upset(fromList(listInput), nsets = 5, sets = rev(colnames(fromList(listInput)))[7:11], cutoff = 10,
      keep.order = TRUE, empty.intersections = "on", nintersects = 5, #group.by = "sets", #order.by = c("freq", "degree"),
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

# Set scale range to be consistent across all custom pHeatMaps
breaksList = seq(-0.2, 16, by = 1)

#------------------------------------------------------------------------------
# Heatmap 3B
# Comparison 1: res_uterus_vs_liver
#------------------------------------------------------------------------------

favourites <- c("ENSMUSG00000062960", "ENSMUSG00000029231", "ENSMUSG00000047414", "ENSMUSG00000006369", "ENSMUSG00000017737",
                "ENSMUSG00000035493", "ENSMUSG00000022150", "ENSMUSG00000039239", "ENSMUSG00000041324", "ENSMUSG00000052187",
                "ENSMUSG00000055609", "ENSMUSG00000052217", "ENSMUSG00000072115", "ENSMUSG00000004328", "ENSMUSG00000027993", 
                "ENSMUSG00000001507" )

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
matrix2                <- as.matrix(mat2.ann)
matrix2.means          <- as.matrix(mat2.ann.means)

ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(length(breaksList))
AnnotCols <- list( tissue=c("liver"="lightgrey", "uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"))

pdf(paste0(Project, "_Fig.3B.allsamples.pdf"), width=7.5,height=12.5, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, breaks = breaksList, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2, col=ScaleCols, cutree_cols=3, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.tissue/res_uterus_vs_liver"))
dev.off()

pdf(paste0(Project, "_Fig.3B.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2.means, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2.means, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.tissue/res_uterus_vs_liver"))
dev.off()




#------------------------------------------------------------------------------
# Heatmap 3C
# Comparison 2: res_uterus_trNK_uterus_ILC1_vs_liver_ILC1
# c("02", "04", "07", "08", "09", "12", "13", "14")
#------------------------------------------------------------------------------

favourites <- c("ENSMUSG00000059256", "ENSMUSG00000022156", "ENSMUSG00000015441", "ENSMUSG00000040284", "ENSMUSG00000021939",
                "ENSMUSG00000038642", "ENSMUSG00000025856", "ENSMUSG00000005413", "ENSMUSG00000021822", "ENSMUSG00000049723",
                "ENSMUSG00000023951", "ENSMUSG00000021594", "ENSMUSG00000027068", "ENSMUSG00000032122", "ENSMUSG00000025464",
                "ENSMUSG00000035385", "ENSMUSG00000035352", "ENSMUSG00000029304", "ENSMUSG00000035373", "ENSMUSG00000009185",
                "ENSMUSG00000031780", "ENSMUSG00000040899", "ENSMUSG00000050232", "ENSMUSG00000030336", "ENSMUSG00000030156",
                "ENSMUSG00000048521", "ENSMUSG00000024610")

res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.df                 <- as.data.frame(res_uterus_trNK_uterus_ILC1_vs_liver_ILC1)
res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.df                 <- res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.df[ favourites, ]
res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.df$ensembl_gene_id <- rownames(res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.df)

rows                       <- match(res_uterus_trNK_uterus_ILC1_vs_liver_ILC1.df$ensembl_gene_id, row.names(rld.group))
mat2                       <- assay(rld.group)[rows,]
df2                        <- as.data.frame(colData(rld.group)[,c("tissue", "collated1")])
colnames(df2)              <- c("tissue", "Comparison")
df2                        <- df2[c("02", "04", "07", "08", "09", "12", "13", "14"),]
df2$Comparison             <- gsub("utrNK_uILC1", "uterus_trNK_ILC1", df2$Comparison)
df2$Comparison             <- gsub("lILC1", "liver_ILC1", df2$Comparison)
df2.means                  <- as.data.frame( list( "Group"=c("uterus_trNK_ILC1", "liver_ILC1"),  "tissue"=c("uterus","liver")) )
rownames(df2.means)        <- c("uterus_trNK_ILC1", "liver_ILC1")
mat2.df                    <- as.data.frame(mat2)
mat2.df$ensembl_gene_id    <- rownames(mat2.df)
mat2.ann                   <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
mat2.ann                   <- mat2.ann[order(sapply(mat2.ann$ensembl_gene_id, function(x) which(x == favourites))), ]
rownames(mat2.ann)         <- mat2.ann$external_gene_name
mat2.ann                   <- within(mat2.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))
mat2.ann.means             <- mat2.ann
mat2.ann.means$uterus_trNK_ILC1 <- rowMeans(subset(mat2.ann, select = c("02", "07", "12", "08","13")), na.rm = TRUE)
mat2.ann.means$liver_ILC1  <- rowMeans(subset(mat2.ann, select = c("04","09","14")),   na.rm = TRUE)
mat2.ann.means             <- mat2.ann.means[ , c("uterus_trNK_ILC1", "liver_ILC1") ]
mat2.ann                   <- mat2.ann[, c("02", "04", "07", "08", "09", "12", "13", "14")]
matrix2                    <- as.matrix(mat2.ann)
matrix2.means              <- as.matrix(mat2.ann.means)

ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(length(breaksList))
AnnotCols       <- list( tissue=c("liver"="lightgrey", "uterus"="darkgrey"),  cell = c("ILC1"="#EA5160", "trNK"="#0099FF"), Comparison=c("uterus_trNK_ILC1"="#0099FF", "liver_ILC1"="#EA5160"))
AnnotCols.means <- list( tissue=c("liver"="lightgrey", "uterus"="darkgrey"),  cell = c("ILC1"="#EA5160", "trNK"="#0099FF"), Group=c("uterus_trNK_ILC1"="#0099FF", "liver_ILC1"="#EA5160"))

pdf(paste0(Project, "_Fig.3C.allsamples.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_trNK_uterus_ILC1_vs_liver_ILC1"))
dev.off()

pdf(paste0(Project, "_Fig.3C.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2.means, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2.means, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols.means, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_trNK_uterus_ILC1_vs_liver_ILC1"))
dev.off()




#------------------------------------------------------------------------------
# Heatmap 3D
# Comparison 3: res_uterus_ILC1_uterus_trNK_vs_uterus_cNK
#  c("08","13",  "02","07","12",  "01","06")
#------------------------------------------------------------------------------

favourites <- c("ENSMUSG00000002985","ENSMUSG00000008540","ENSMUSG00000060063","ENSMUSG00000004207","ENSMUSG00000000817",
                "ENSMUSG00000024401","ENSMUSG00000050395","ENSMUSG00000025498","ENSMUSG00000018774","ENSMUSG00000031447",
                "ENSMUSG00000016534","ENSMUSG00000052688","ENSMUSG00000035373","ENSMUSG00000035385","ENSMUSG00000044827",
                "ENSMUSG00000050240","ENSMUSG00000025746","ENSMUSG00000028470","ENSMUSG00000030147","ENSMUSG00000045087",
                "ENSMUSG00000079620")

res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.df                 <- as.data.frame(res_uterus_ILC1_uterus_trNK_vs_uterus_cNK)
res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.df                 <- res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.df[ favourites, ]
res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.df$ensembl_gene_id <- rownames(res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.df)

rows                       <- match(res_uterus_ILC1_uterus_trNK_vs_uterus_cNK.df$ensembl_gene_id, row.names(rld.group))
mat2                       <- assay(rld.group)[rows,]
df2                        <- as.data.frame(colData(rld.group)[,c("tissue", "collated2")])
colnames(df2)              <- c("tissue", "Comparison")
df2                        <- df2[c("08","13",  "02","07","12",  "01","06"), ]

df2$Comparison             <- gsub("uILC1_utrNK", "uterus_ILC1_trNK", df2$Comparison)
df2$Comparison             <- gsub("ucNK", "uterus_cNK", df2$Comparison)

df2.means                  <- as.data.frame( list( "Comparison"=c("uterus_ILC1_trNK", "uterus_cNK"), "tissue"=c("uterus","uterus")) )
rownames(df2.means)        <- c("uterus_ILC1_trNK", "uterus_cNK")
mat2.df                    <- as.data.frame(mat2)
mat2.df$ensembl_gene_id    <- rownames(mat2.df)
mat2.ann                   <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
mat2.ann                   <- mat2.ann[order(sapply(mat2.ann$ensembl_gene_id, function(x) which(x == favourites))), ]
rownames(mat2.ann)         <- mat2.ann$external_gene_name
mat2.ann                   <- within(mat2.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))
mat2.ann.means             <- mat2.ann
mat2.ann.means$uterus_ILC1_trNK <- rowMeans(subset(mat2.ann, select = c("08","13", "02", "07", "12")),        na.rm = TRUE)
mat2.ann.means$uterus_cNK       <- rowMeans(subset(mat2.ann, select = c("01", "06")),   na.rm = TRUE)
mat2.ann.means             <- mat2.ann.means[ , c("uterus_ILC1_trNK", "uterus_cNK") ]
mat2.ann                   <- mat2.ann[, c("08","13",  "02","07","12",  "01","06")]
matrix2                    <- as.matrix(mat2.ann)
matrix2.means              <- as.matrix(mat2.ann.means)

ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(length(breaksList))
AnnotCols       <- list( tissue=c("uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"), Comparison=c("uterus_ILC1_trNK"="#0099FF", "uterus_cNK"="#FFCC00"))
AnnotCols.means <- list( tissue=c("uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"), Comparison=c("uterus_ILC1_trNK"="#0099FF", "uterus_cNK"="#FFCC00"))

pdf(paste0(Project, "_Fig.3D.allsamples.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, breaks = breaksList, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_ILC1_uterus_trNK_vs_uterus_cNK"))
dev.off()

pdf(paste0(Project, "_Fig.3D.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2.means, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2.means, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols.means, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_ILC1_uterus_trNK_vs_uterus_cNK"))
dev.off()



#------------------------------------------------------------------------------
# Heatmap 3E
# Comparison 4: res_uterus_cNK_vs_uterus_trNK
#  c("01", "06",  "02", "07", "12")
#------------------------------------------------------------------------------

favourites <- c("ENSMUSG00000039217","ENSMUSG00000022965","ENSMUSG00000055170","ENSMUSG00000026070","ENSMUSG00000026068",
                "ENSMUSG00000000440","ENSMUSG00000020609","ENSMUSG00000015568","ENSMUSG00000025044","ENSMUSG00000002602",
                "ENSMUSG00000014361","ENSMUSG00000040249","ENSMUSG00000007613","ENSMUSG00000035493","ENSMUSG00000032902",
                "ENSMUSG00000027398","ENSMUSG00000070390","ENSMUSG00000032691","ENSMUSG00000026072","ENSMUSG00000031780",
                "ENSMUSG00000009185","ENSMUSG00000042284","ENSMUSG00000015533")

res_uterus_cNK_vs_uterus_trNK.df                 <- as.data.frame(res_uterus_cNK_vs_uterus_trNK)
res_uterus_cNK_vs_uterus_trNK.df                 <- res_uterus_cNK_vs_uterus_trNK.df[ favourites, ]
res_uterus_cNK_vs_uterus_trNK.df$ensembl_gene_id <- rownames(res_uterus_cNK_vs_uterus_trNK.df)

rows                       <- match(res_uterus_cNK_vs_uterus_trNK.df$ensembl_gene_id, row.names(rld.group))
mat2                       <- assay(rld.group)[rows,]
df2                        <- as.data.frame(colData(rld.group)[,c("cell", "tissue")])
df2                        <- df2[c("01", "06",  "02", "07", "12"), ]
df2.means                  <- as.data.frame( list( "Comparison"=c("uterus_cNK", "uterus_trNK"), "tissue"=c("uterus","uterus")) )
rownames(df2.means)        <- c("uterus_cNK", "uterus_trNK")
mat2.df                    <- as.data.frame(mat2)
mat2.df$ensembl_gene_id    <- rownames(mat2.df)
mat2.ann                   <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
mat2.ann                   <- mat2.ann[order(sapply(mat2.ann$ensembl_gene_id, function(x) which(x == favourites))), ]
rownames(mat2.ann)         <- mat2.ann$external_gene_name
mat2.ann                   <- within(mat2.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))
mat2.ann.means             <- mat2.ann
mat2.ann.means$uterus_cNK <- rowMeans(subset(mat2.ann, select = c("01", "06")),        na.rm = TRUE)
mat2.ann.means$uterus_trNK <- rowMeans(subset(mat2.ann, select = c("02", "07", "12")), na.rm = TRUE)
mat2.ann.means             <- mat2.ann.means[ , c("uterus_cNK", "uterus_trNK") ]
mat2.ann                   <- mat2.ann[, c("01", "06",  "02", "07", "12")]
matrix2                    <- as.matrix(mat2.ann)
matrix2.means              <- as.matrix(mat2.ann.means)

ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(length(breaksList))
AnnotCols       <- list( tissue=c("uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "trNK"="#0099FF"))
AnnotCols.means <- list( tissue=c("uterus"="darkgrey"),  cell = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF"), Comparison=c("uterus_trNK"="#0099FF","uterus_cNK"="#FFCC00"))

pdf(paste0(Project, "_Fig.3E.allsamples.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_cNK_vs_uterus_trNK"))
dev.off()

pdf(paste0(Project, "_Fig.3E.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2.means, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2.means, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols.means, cellwidth = 20, fontsize=18, cellheight = 20, main=paste0(Project, "\nrld.group/res_uterus_cNK_vs_uterus_trNK"))
dev.off()



#------------------------------------------------------------------------------
# Heatmap 3F
# Comparison 5: res_uterus_ILC1_vs_uterus_trNK
#  c("08","13",  "02", "07", "12")
#------------------------------------------------------------------------------

favourites <- c("ENSMUSG00000036594","ENSMUSG00000073421","ENSMUSG00000060586","ENSMUSG00000079547","ENSMUSG00000024610",
                "ENSMUSG00000026012","ENSMUSG00000026117","ENSMUSG00000022657","ENSMUSG00000005696","ENSMUSG00000073494",
                "ENSMUSG00000026573","ENSMUSG00000032021","ENSMUSG00000032035","ENSMUSG00000020689","ENSMUSG00000020589",
                "ENSMUSG00000024940","ENSMUSG00000028965","ENSMUSG00000044338","ENSMUSG00000032915","ENSMUSG00000020681",
                "ENSMUSG00000022528")

res_uterus_ILC1_vs_uterus_trNK.df                 <- as.data.frame(res_uterus_ILC1_vs_uterus_trNK)
res_uterus_ILC1_vs_uterus_trNK.df                 <- res_uterus_ILC1_vs_uterus_trNK.df[ favourites, ]
res_uterus_ILC1_vs_uterus_trNK.df$ensembl_gene_id <- rownames(res_uterus_ILC1_vs_uterus_trNK.df)

rows                       <- match(res_uterus_ILC1_vs_uterus_trNK.df$ensembl_gene_id, row.names(rld.group))
mat2                       <- assay(rld.group)[rows,]
df2                        <- as.data.frame(colData(rld.group)[,c("cell", "tissue")])
df2                        <- df2[c("08","13",  "02", "07", "12"), ]
df2.means                  <- as.data.frame( list( "Comparison"=c("uterus_ILC1", "uterus_trNK"), "tissue"=c("uterus","uterus")) )
rownames(df2.means)        <- c("uterus_ILC1", "uterus_trNK")
mat2.df                    <- as.data.frame(mat2)
mat2.df$ensembl_gene_id    <- rownames(mat2.df)
mat2.ann                   <- merge(mat2.df, ensEMBL2id, by="ensembl_gene_id")
mat2.ann                   <- mat2.ann[order(sapply(mat2.ann$ensembl_gene_id, function(x) which(x == favourites))), ]
rownames(mat2.ann)         <- mat2.ann$external_gene_name
mat2.ann                   <- within(mat2.ann, rm("ensembl_gene_id", "external_gene_name", "description", "entrezgene", "chromosome_name", "gene_biotype", "Length"))
mat2.ann.means             <- mat2.ann
mat2.ann.means$uterus_ILC1 <- rowMeans(subset(mat2.ann, select = c("08", "13")),       na.rm = TRUE)
mat2.ann.means$uterus_trNK <- rowMeans(subset(mat2.ann, select = c("02", "07", "12")), na.rm = TRUE)
mat2.ann.means             <- mat2.ann.means[ , c("uterus_ILC1", "uterus_trNK") ]
mat2.ann                   <- mat2.ann[, c("08","13",  "02", "07", "12")]
matrix2                    <- as.matrix(mat2.ann)
matrix2.means              <- as.matrix(mat2.ann.means)

ScaleCols <- colorRampPalette(colors = c("yellow","red","black"))(length(breaksList))
AnnotCols       <- list( tissue=c("uterus"="darkgrey"),  cell = c("ILC1"="#EA5160", "trNK"="#0099FF"))
AnnotCols.means <- list( tissue=c("uterus"="darkgrey"),  cell = c("ILC1"="#EA5160", "trNK"="#0099FF"), Comparison=c("uterus_trNK"="#0099FF","uterus_ILC1"="#EA5160"))

pdf(paste0(Project, "_Fig.3E.allsamples.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2, cutree_cols=3, col=ScaleCols, annotation_colors=AnnotCols, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_ILC1_vs_uterus_trNK"))
dev.off()

pdf(paste0(Project, "_Fig.3E.pdf"), width=7.5,height=10, onefile=FALSE)
par(bg=NA)
pheatmap(matrix2.means, breaks = breaksList, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE, cluster_cols=TRUE, annotation_col=df2.means, cutree_cols=2, col=ScaleCols, annotation_colors=AnnotCols.means, cellwidth = 20, cellheight = 20, fontsize=18, main=paste0(Project, "\nrld.group/res_uterus_ILC1_vs_uterus_trNK"))
dev.off()

#
# pHeatMap Figure 2A
#

logFoldChanceCutOff <- 7.5
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
#mat3.ann                <- within(mat3.ann, rm("ensembl_gene_id","external_gene_name","description","entrezgene","chromosome_name","gene_biotype","Length"))
mat3.ann                <- within(mat3.ann, rm("ensembl_gene_id","external_gene_name","description","entrezgene","chromosome_name","gene_biotype"))
matrix3                 <- as.matrix(mat3.ann)

colors<-colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(255)

pdf(paste0(Project, "_Fig.2A.l2fc.", logFoldChanceCutOff, ".padj.", significance,  ".pdf"), width=7,height=12, onefile=FALSE)
par(bg=NA)
pheatmap(matrix3, cluster_rows=TRUE, show_rownames=FALSE, show_colnames=FALSE, cluster_cols=TRUE, annotation_col=df3, 
         annotation_colors=list(tissue=c(liver="lightgrey",uterus="darkgrey"), cell=c(cNK="orange", ILC1="firebrick3", trNK="steelblue3")),
         cutree_cols=6, border_color="white",
         col=colors, main=paste(Project, " pHeatMap \nrld.groups [log2FC=", logFoldChanceCutOff, " padj=", significance, "]",  sep=""))
dev.off()








message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

customPCA <- function(sampleTBL, RLD, TOPNUM, model) {
  
  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, condition=sampleTBL$group, Tissue=sampleTBL$tissue, Cell=sampleTBL$cell)
   
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=Cell, fill=Cell, shape=(factor(scores$Tissue)), label=scores$sampleName) ) +
             geom_point(size = 3, alpha=0.75 ) + 
             geom_text_repel(aes(label=sampleName), show.legend = FALSE) +
             xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
             ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
             scale_shape_manual(name="Tissue", values = c(17, 16)) + 
             scale_fill_manual(name="Cell Type",  values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +
             scale_colour_manual(name="Cell Type", values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +    
             theme(text = element_text(size=elementTextSize)) 
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=Cell, fill=factor(scores$Cell), shape=(factor(scores$Tissue)), label=scores$sampleName) ) +
                #geom_encircle(alpha = 0.2, show.legend = FALSE) +
                geom_point(size = 6, alpha=0.75 ) + 
                xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
                ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
                scale_shape_manual(name="Tissue", values = c(17, 16)) + 
                scale_fill_manual(name="Cell Type",  values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +
                scale_colour_manual(name="Cell Type", values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +    
                theme(text = element_text(size=elementTextSize)) 
  
  dd.row          <- as.dendrogram(hclust(dist(t( RLD[select, ]))))
  ddata           <- dendro_data(dd.row)
  labs            <- label(ddata)
  labs$sampleName <- labs$label
  annotations     <- merge(labs, sampleTBL, by="sampleName")
  
  plt.dendro <- ggplot(segment(ddata)) +
                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                geom_text(data=annotations,aes(label=paste0(cell, "-", tissue), x=x, y=-0.25, angle=-90, colour=annotations$cell), hjust=0) +
                scale_colour_manual(name="Cell Type", values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +   
                scale_y_continuous(breaks=seq(0, max(ddata$segments$y), 50)) +
                expand_limits(y=c( -1*(max(ddata$segments$y)/2)), max(ddata$segments$y) ) + 
                ylab("distance") + xlab("") +
                ggtitle(paste0("Dendrogram Top ", TOPNUM, " MV")) +
                theme(text = element_text(size=elementTextSize), legend.position="none",
                      axis.text.x = element_blank(),  axis.ticks.x = element_blank() )
  
  return(list(plt.pca, plt.pca.nl, plt.dendro) )
}


pca.plt.rld.250            <- customPCA(sampleTable, assay(rld.group), 250, "rld")
pca.plt.rld.250[[1]]
pca.plt.rld.250[[2]]

pca.plt.rld.500            <- customPCA(sampleTable, assay(rld.group), 500, "rld")
pca.plt.rld.1000           <- customPCA(sampleTable, assay(rld.group), 1000, "rld")
pca.plt.rld.47000          <- customPCA(sampleTable, assay(rld.group), 47000, "rld")

pdf(paste(Project, "_Fig.PCA.pdf", sep=""),width=10,height=7)
par(bg=NA)
plot_grid(pca.plt.rld.250[[2]], pca.plt.rld.500[[2]], pca.plt.rld.1000[[2]], pca.plt.rld.47000[[2]], nrow=2, ncol=2)
dev.off()

pdf(paste(Project, "_Fig.PCA_250MV.pdf", sep=""),width=10,height=7)
par(bg=NA)
pca.plt.rld.250[[2]]
dev.off()

pdf(paste(Project, "_Fig.hclust.pdf", sep=""),width=10,height=7)
par(bg=NA)
plot_grid(pca.plt.rld.250[[3]], pca.plt.rld.500[[3]], pca.plt.rld.1000[[3]], pca.plt.rld.47000[[3]], nrow=2, ncol=2)
dev.off()


library("RColorBrewer")
sampleDists                <- dist(t(assay(vsd.group)))
sampleDistMatrix           <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.group$tissue, vsd.group$cell, sep="-")
colnames(sampleDistMatrix) <- NULL
colors                     <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste(Project, "_Fig.SampleSimMatrix.pdf", sep=""),width=8,height=7)
par(bg=NA)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()




rv     = rowVars(assay(rld.group))
select = order(rv, decreasing = TRUE)[seq_len(min(47000, length(rv)))]
pca    = prcomp(t(assay(rld.group)[select, ]))
fac    = factor(apply(as.data.frame(colData(rld.group)[, "group", drop = FALSE]), 1, paste, collapse = " : "))
                        
pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(sampleFiles=sampleTable$fileName, pca$x, sampleCell=sampleTable$cell, sampleName=sampleTable$sampleName, sampleTissue=sampleTable$tissue)

pdf(paste(Project, "_Fig.2B.pdf", sep=""),width=10,height=7)
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

pdf(paste(Project, "_Fig.2B.labels.pdf", sep=""),width=10,height=7)
par(bg=NA)
ggplot(data=scores, 
       aes(x = scores$PC1, y = scores$PC2, col = factor(scores$sampleCell), shape=(factor(scores$sampleTissue)), label=scores$sampleFiles )) +
  geom_point(size = 5 ) + 
  geom_text_repel(aes(label=scores$sampleName), box.padding = unit(1.0, "lines")) +
  xlab(pc1lab) + ylab(pc2lab) + ggtitle(paste(Project, " PCA (All Genes)", sep="")) +
  scale_shape_manual(name="Tissue", values = c(17, 16)) + 
  scale_fill_manual(name="Cell Type",  values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +
  scale_colour_manual(name="Cell Type", values = c("cNK"="#FFCC00", "ILC1"="#EA5160", "trNK"="#0099FF")) +
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

pdf(paste(Project, "_SuppFig.S2A.pdf", sep=""),width=10,height=7)
par(bg=NA)
pca.1.25.plot
dev.off()

pca.2         <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <-  pca.2[c(1:25),]
pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pdf(paste(Project, "_SuppFig.S2B.pdf", sep=""),width=10,height=7)
par(bg=NA)
pca.2.25.plot
dev.off()






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
  
  pdf(paste(Project, "_Figure.", SetName, ".pdf", sep=""),width=3.5,height=myHeight, onefile=FALSE)
  par(bg=NA)
  print({ plot })
  dev.off()
  
  return(plot)
}


# Figure 2C p1    2Cp1
selected.2Cp1    <- c("ENSMUSG00000015533", "ENSMUSG00000030325", "ENSMUSG00000042284", "ENSMUSG00000062524")
plt.2Cp1 <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.2Cp1,"2Cp1")
# Figure 2C p2    2Cp2
selected.2Cp2    <- c("ENSMUSG00000001444", "ENSMUSG00000032446", "ENSMUSG00000055170")
plt.2Cp2 <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.2Cp2,"2Cp2")
# Figure 2D       2D
selected.2D    <- c("ENSMUSG00000004730", "ENSMUSG00000014361","ENSMUSG00000022901", "ENSMUSG00000026395")
plt.2D   <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.2D,"2D")


# Figure 4Ap1 
selected.4A    <- c("ENSMUSG00000048521", "ENSMUSG00000030173")
plt.4A <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.4Ap1,"4A")

# Figure 4Fp1 
selected.4Fp1    <- c("ENSMUSG00000030173", "ENSMUSG00000079852", "ENSMUSG00000079853", "ENSMUSG00000089727")
plt.4Fp1 <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.4Fp1,"4Fp1")

# Figure 4Fp2 
selected.4Fp2    <- c("ENSMUSG00000033024", "ENSMUSG00000067599", "ENSMUSG00000005947", "ENSMUSG00000026581", "ENSMUSG00000031004")
plt.4Fp2 <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.4Fp2,"4Fp2")

# Figure 4Fp3 
selected.4Fp3    <- c("ENSMUSG00000030149", "ENSMUSG00000030167","ENSMUSG00000030156", "ENSMUSG00000030336")
plt.4Fp3 <- makeGeneSetExpressionPlot(dds.group.FPKM,Project,selected.4Fp3,"4Fp3")





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


CairoPDF(file = paste0(Project, "_Fig.1B.pdf"), width = 12, height = 6, onefile = TRUE, family = "Arial" )
par(bg=NA)
all.points.plot
dev.off()

CairoPDF(file = paste0(Project, "_SuppFig.S1.pdf"), width = 12, height = 6, onefile = TRUE, family = "Arial" )
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

CairoPDF(file = paste0(Project, "_Fig.S1.inset.pdf"), width = 4, height = 6, onefile = TRUE, family = "Arial" )
par(bg=NA)
inset.df.plot
dev.off()

CairoPDF(file = paste0(Project, "_SuppFig.S1.inset.pdf"), width = 4, height = 6, onefile = TRUE, family = "Arial" )
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
Immune_cells_expr_matrix                <- read.csv( paste0(Base.dir, "/DataTables/srep40508-s1.csv"))
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
           theme(text = element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")

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


#pdf(paste0(Project, "-ImmuCC_DeconRNA-Seq_Grid.pdf"),width=20,height=20, onefile=FALSE)
#par(bg=NA)
#plot_grid(plot.Mast.Cells, plot.Neutrophil.Cells, plot.Eosinophil.Cells, plot.B.Cells.Memory, plot.B.Cells.Naive,
#          plot.Plasma.Cells, plot.T.Cells.CD8.Actived, plot.T.Cells.CD8.Naive, plot.T.Cells.CD8.Memory,plot.M0.Macrophage,
#          plot.M1.Macrophage, plot.M2.Macrophage, plot.Treg.Cells, plot.T.Cells.CD4.Memory, plot.T.Cells.CD4.Naive,
#          plot.T.Cells.CD4.Follicular, plot.Th1.Cells, plot.Th17.Cells, plot.Th2.Cells, plot.Monocyte,
#          plot.GammaDelta.T.Cells, plot.NK.Resting, plot.NK.Actived, plot.DC.Actived, plot.DC.Immature,
#          ncol = 5, nrow = 5)
#dev.off()


head(Decon_results_df1)
head(Decon_results_df2)
head(Decon_results_df3)

Decon_results_df <- merge(Decon_results_df1, Decon_results_df2,by="row.names")
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL

Decon_results_df <- merge(Decon_results_df,  Decon_results_df3,by="row.names",all.x=TRUE)
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL


Decon_results_df.plot <- Decon_results_df
Decon_results_df.plot <- within(Decon_results_df, rm("sampleName", "cell.x", "tissue.x", "group.x", "cell.y", "tissue.y", "group.y"))
Decon_results_df.plot <- melt(Decon_results_df.plot)
head(Decon_results_df.plot)

pdf(paste0(Project, "_Fig.S4A.pdf"),width=12.5,height=15, onefile=FALSE)
par(bg=NA)
ggplot(Decon_results_df.plot, aes(x = group, y = value)) +
  geom_boxplot(width = 0.1, alpha=0.5) +
  geom_point(position=position_jitter(w=0.05,h=0), alpha=0.75, colour="red") +
  labs(y = "Cell type fraction", x = "") +
  facet_wrap( ~ variable , ncol=5) +
  theme_minimal() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle = 90),
        legend.position="none",
        strip.text.x = element_text(face='bold'),
        strip.background = element_rect(colour="white", fill='white'))
dev.off()


Decon_results_df.means      <- aggregate(Decon_results_df[, c(1:10,15:28)], list(Decon_results_df$group), mean)
Decon_results_df.means.melt <- melt(Decon_results_df.means)

Decon_results_df.vars      <- aggregate(Decon_results_df[, c(1:10,15:28)], list(Decon_results_df$group), var)
Decon_results_df.vars.melt <- melt(Decon_results_df.vars)

Decon_results_df.means.melt$vars      <- Decon_results_df.vars.melt$value

Decon_results_df.means.melt$normvars  <- Decon_results_df.means.melt$vars/max(Decon_results_df.means.melt$vars)

colnames(Decon_results_df.means.melt) <- c("Group", "Immune.Cell", "mean", "variance", "normalisedvariance")


deconRNASeq.summary.table <- Decon_results_df.means.melt
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


pdf(paste0(Project, "_Fig.S4B.pdf"),width=20,height=20, onefile=FALSE)
par(bg=NA)
ggplot(deconRNASeq.summary.table, aes(y = Immune.Cell,x = Group)) +
  geom_point(aes(size = mean, colour = Group, alpha=(1-normalisedvariance))) + 
  xlab("") + ylab("") +
  guides(alpha=guide_legend(title="1-(Normalised Variance)")) +
  guides(size=guide_legend(title="Mean Proportion")) +
  guides(colour=guide_legend(title="Group")) +  
  guides(colour=FALSE) +
  theme_bw()
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ Analysis Complete")
message("+-------------------------------------------------------------------------------")


#------------------------------------------------------------------------------
# FIN
#------------------------------------------------------------------------------

