# 2018-Filipovic-Colucci [TITLE TO BE DETERMINED]

**Iva Filipovic<sup>1,2,3</sup>, Russell S. Hamilton<sup>2,3</sup>, [OTHER AUTHORS], Francesco Colucci<sup>1,3,§</sup>**

<sup>1</sup> Department of Obstetrics and Gynaecology, University of Cambridge, University of Cambridge School of Clinical Medicine, NIHR Cambridge Biomedical Research Centre <br>
<sup>2</sup> Department of Physiology, Development and Neuroscience, University of Cambridge<br>
<sup>3</sup> Centre for Trophoblast Research, University of Cambridge<br>
<sup>§</sup> Corresponding author: fc287@medschl.cam.ac.uk <br>

### Abstract ###
To be added

### Data Processing ###
Data were aligned to GRCm38 mouse genome (Ensembl Release 84) with TopHat2 (v2.1.1, using bowtie2 v2.2.9) with a double map strategy. Alignments and QC were processed using custom ClusterFlow (v0.5dev) pipelines and assessed using MultiQC (0.9.dev0). Gene quantification was determined with HTSeq-Counts (v0.6.1p1). Additional quality control was performed with feature counts (v 1.5.0-p2), qualimap (v2.2) and preseq (v2.0.0). Differential gene expression was performed with DESeq2 package (v1.18.1, R v3.4.2) and with the same package read counts were normalised on the estimated size factors.

Resource       | URL
-------------- | --------------
GRCm38         | [Link](http://mar2016.archive.ensembl.org/index.html)
FastQC         | [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Trim_galore    | [Link](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
TopHat2        | [DOI](http://dx.doi.org/10.1186/gb-2013-14-4-r36)
ClusterFlow    | [DOI](http://dx.doi.org/10.12688/f1000research.10335.2)
MultiQC        | [DOI](http://dx.doi.org/10.1093/bioinformatics/btw354)
HTSeq-counts   | [DOI](http://dx.doi.org/10.1093/bioinformatics/btu638)
Feature_counts | [DOI](http://dx.doi.org/10.1093/bioinformatics/btt656)
Qualimap       | [DOI](https://doi.org/10.1093/bioinformatics/bts503)

A custom module for TopHat2 double map is provided in this repository, and can be run, by copying it into the modules directory of a ClusterFlow installation. With HTSeq-Counts gene count tables figures in the table below can be reproduced with the R script provided in this repository.

### Script to reproduce paper figures ###

The provided R script assumes the script is placed in a directory containing a subdirectory (called HTSeq_Counts) with all the htseq-counts files (one per sample). The script can be run interactively in R-studio or as a batch using Rscript. Note that some of the figures in the manuscript have had some label positions moved manually to prevent overlaps.

Figure    | Description | Output Filename
--------- | ----------- | ------------------------
Figure 1  | Some Plot   | 2018-Filipovic-Colucci_Fig1.pdf


_Additional Data_

Description                    | Output Filename
------------------------------ | ------------------------
DEG Results: Uterine_vs_Liver  | 2018-Filipovic-Colucci_DESeq2_DEGs_Uterine_vs_Liver.csv

### Sample Table ###

sample Name  | cell | tissue |
------------ | ---- | ------ |   
SLX-9343.01  | cNK  | uterus  
SLX-9343.02  | trNK | uterus
SLX-9343.04  | ILC1 | liver  
SLX-9343.05  |  cNK | liver   
SLX-9343.06  | cNK  | uterus  
SLX-9343.07  | trNK | uterus
SLX-9343.08  | ILC1 | uterus
SLX-9343.09  | ILC1 | liver
SLX-9343.10  | cNK  | liver   
SLX-9343.12  | trNK | uterus
SLX-9343.13  | ILC1 | uterus
SLX-9343.14  | ILC1 | liver  
SLX-9343.15  | cNK  | liver   

### Links ###

Description   | URL
------------- | ----------
Publication   | [Journal](http://) and [DOI](http://) <br> (To be updated on publication)
Raw Data      | ArrayExpress EMBL-EBI <br>Data to be released on publication<br> [E-MTAB-XXXX](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5803)
Colucci Group | [Colucci group website](http://moffettcoluccilab.org/francesco-colucci/)

### Contact ###

Contact rsh46 -at- cam.ac.uk for bioinformatics related queries
