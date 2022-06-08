# Project Description

Following the methods of Wang et al in _A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data_ we sought out to recreate the analysis from this paper, which was to determine the concordance of microarray and RNA-seq technologies. Data was provided as a part of BF528.

# Contributors

David Lenci, Daniel Gealow, and Nikita Tomar

# Repository Contents

Qsub scripts directory contain scripts used to perform processing prior to analysis, like running featureCounts that were submitted as SCC jobs. R-scripts and python scripts contain code used to perform analysis on our generated data from the qsub scripts. This includes DESeq2, limma, and concordance calculations. Within the R-scripts both DE_analysis.R and combine_csv.R are ran using the associated R-markdown file programmer_report.RMD, which containes additional instruction/clarifications. 
