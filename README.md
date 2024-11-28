# Cut-RunAnalysis


This python script can be used to analyze Cut&Run or ChIP Seq data for changes in histone marks or transcription factor binding on the X chromosome before and after triggering Xist expression. It reads in Bamfiles, and generates an analysis using genome annotation in GTF format and HTSeq read count data.

Before running the analysis, it is recommended that the generated Bamfiles are filtered using lists of problematic regions that will result in artifactual peaks. Such a list is for example provided in Nordin et al., 2023 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03027-3#MOESM2). 

Output of Xchipanalysis are coverage Bedgraph files for each sample and a dictionary containing statistics of the experiment.


## Usage
To perform analysis, the following parameters have to be specified starting on line 814:

CHIPSEQ_ANALYSIS: True/False 

LOAD_CURRENT_ANALYSIS: True to load in an analysis that has been generated previously and is stored in the DATA_FILE_PATH; otherwise False.

ANNOTATION_FILE_PATH: Location where annotation files are stored. In the download folder, users can find the files built for the mouse genome according to the assembly of Dec. 2011 (GCRm38/mm10).

DATA_FILE_PATH: Location where the analysis should be stored.

BAM_FILE_PATH: Location of bam files that are read in. 

HASAPPY_ANNOTATION_PICKLE_FILE_PATH: location of genome annotation pickl file in HaSAPPY format.

Experiments are added by using 
```Python
analysis.addExperiment("Name_sample", BAM_FILE_PATH+'input.bam', aQualValue=10, frgmtLength=200, maxReadsPerPos=5 )
```

Ater running the script, heatmaps can be generated for various aspects of the dataset by calling the following functions:
```Python
plotFeaturesHeatmaps()
```
generates heatmaps and cumulative coverage plots for the following features: ESs H3K27Me3 peaks, Cpg Islands, ESC Enhancers, ESC Ezh2 binding sites (strong and moderate), TS DNase I Hypersensitivity Sites (DHS), TS Formaldehyde-Assisted Isolation of Regulatory Elements (FAIRE) sites, TS CTCF sites and TS Polymerase II binding sites.
```Python
plotAffectedGeneGroups()
```
generates heatmaps and cumulative coverage plots for two classes of genes on the X chromosomes; "top genes", which are genes whose expression changed upon the addition of doxycylin and the upregulation of Xist RNA in ESCS, and "flop genes", genes whose expression was not greatly affected by Xist RNA expression upon doxycyling treatment, according to the RNA sequencing perfomed in Monfort et al., 2024.
```Python
plotExpressionLevelGroups()
```
generates heatmaps and cumulative coverage plots for three classes of genes on the X chromosome; genes which were highly, medium and low expressed in ESCs, according to the RNA sequencing perfomed in Monfort et al., 2024.
```Python
analyseGeneGroups()
```
generates cumulative profiles separately for X chromosomal and autosomal genes, and plots them over the gene body and around the transcriptional start site (TSS).