################################################################
### Quick script to filter and reorganize count matrix files ###
################################################################
# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# Readme # ---------------------
# This script will load the count matrix and filter it for low abundant gene counts
# and then selects only the samples form control and high exposure.
# Finally both count matrices are brought to the common set of genes and exported to
# csv. This file formating is needed for Reactome.
##########

filterFUN = function(mtx,anno){
  # select IDs for columns 
  ID = rownames(anno[grep("[Cc]ontrol",anno$Condition),])
  HE = rownames(anno[grep("[Hh]igh",anno$Condition),])
  # filter low read counts
  keep <- names(which(rowSums(mtx) >= ncol(mtx)+2))
  mtx[keep,c(ID,HE)]
}

# import 6PTU
message("Loading data ...")
mtx = read.table("../ArrayExpress/E-MTAB-9054/6PTU_CountMatrix.txt", row.names = 1)
anno = read.csv2("../ArrayExpress/E-MTAB-9054/6PTU_coldata.csv", row.names = 1)
PTU = filterFUN(mtx,anno)

# import T3
mtx1 = read.table("../ArrayExpress/E-MTAB-9056/T3_CountMatrix.txt", row.names = 1)
anno1 = read.csv2("../ArrayExpress/E-MTAB-9056/T3_coldata.csv", row.names = 1)
T3 = filterFUN(mtx1,anno1)

# reduce tmp and tmp1 to common set of genes
int = intersect(row.names(PTU),row.names(T3)) #23989 in common set 
# Export matrix as csv
message("Saving filtered matrix to csv file ...")
write.csv2(PTU[int,],"../ArrayExpress/E-MTAB-9054/6PTU_CountMatrix_red.csv")
write.csv2(T3[int,],"../ArrayExpress/E-MTAB-9056/T3_CountMatrix_red.csv")
message("Done!\n")