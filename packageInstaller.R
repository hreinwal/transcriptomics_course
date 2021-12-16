### Installing required R packages for course ###

# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# README # ---------------------------------------------------------
# Small script to quickly equip the student's R environment with the
# useful and required packages for the ecotox transcriptomics course.
####################################################################

## Bioconductor source packages ##
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BioconPackage = c("DESeq2","IHW","ComplexHeatmap","ReactomePA",
                  "apeglm","qvalue","pheatmap","org.Dr.eg.db",
                  "EnhancedVolcano","pcaMethods","AnnotationDbi",
                  "clusterProfiler","biomaRt","vsn")
BiocManager::install(BioconPackage, ask = F)

## CRAN source packages ##
CranP = c("RColorBrewer","viridis","dplyr","ggplot2","gtools","vegan","eulerr",
          "ggpubr","hexbin","Rtsne","devtools","tidyr","tidyverse","venn")
for(i in CranP) install.packages(i)
rm(i)

### END OF SCRIPT ###
