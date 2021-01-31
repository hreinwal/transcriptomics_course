########################
###   ArrayExpress   ###
########################

# Author:  Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# README # ---------------------------------
# This is a training script to learn the basics of the ArrayExpress R package.
# https://www.bioconductor.org/packages/release/bioc/html/ArrayExpress.html
# ArrayExpress is a large public database for experimental RNAseq data.
# It stores Functional Genomics Data data from high-throughput functional genomics experiments,
# and provides these data for reuse to the research community.
# https://www.ebi.ac.uk/arrayexpress/
# It's entries are linked to the European Nucleotide Archive (ENA) where the raw fastq files are stored
# https://www.ebi.ac.uk/ena/browser/home

# We can use the ArrayExpress package to easily download raw and processed data as well as 
# meta information about the samples via their accession number. With the package this can be
# quickly done in R.
##########

### Install ArrayExpress ###
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("ArrayExpress")

### Load the package ###
require(ArrayExpress)
browseVignettes("ArrayExpress")

### Prepare environment ###
setwd("../") #go one dir up (out of the git repo)
home = getwd()

#############################################
### Download data from ArrayExpress via R ###
#############################################

# In order to download data from ArrayExpress you need to know the accession number of the 
# database entries which you would like to download. 
# You can search the database online through your browser: https://www.ebi.ac.uk/arrayexpress/
# or directly in R with ?queryAE()

# Now we are searching AE for danio rerio with the keyword "thyroid"
sets = queryAE(keywords = "thyroid", species = "danio+rerio")

# Q.1) Which ArrayExpress (AE) entries provide data that was published in a scientific journal?
# Q.2) Which substances were tested in those experiments?
# Q.3) How would you search for entries linked to cancer research in zebrafish? How many entries can you find?
#      Try to retrieve cancer related entries for humans. How many can you find there?


### Import an ArrayExpress dataset into R ### ----------------------------
# Once  you  know  which  identifier  you  wish  to  retrieve,  you can use
# ?ArrayExpress() to load the data into R
rawset = ArrayExpress(accession = "E-MTAB-9056") 
# This bugs! Tells me there is no raw data although there is ...


### use getAE to store Arrayexpress datasets locally ### ------------------

# To organize the output a little bit nice I build these functions.
# AEDownload will generate a folder named after the !single! accession number provided and stores
# the files in there. It will further remove any zipped files.
# If import = TRUE the function will generate an R list object named after the accession containing
# the idf, sdrf,, and processed matrix files.

AEDownload <- function(accession, type = "processed", out = getwd(), import = T) {
  if(length(accession)>1){stop("length(accession) > 1. Please provide only a single accession number. ")}
  # create output dir
  dir.create(paste0(out,"/ArrayExpress"), showWarnings = F)
  setwd(paste0(out,"/ArrayExpress"))
  
  # run getAE() in loop for accession numbers
  i = accession
  dir.create(i, showWarnings = F)
  AE <- getAE(i, type = type, path = i, extract = T)
  zip <- list.files(i, full.names = T)[grep(".zip",list.files(i))]
  if(file.exists(zip[1])) file.remove(zip)
  
  # Import data in R list object
  if(import == T) {
    ls <- list()
    sdrf <- list.files(i, full.names = T)[grep(".sdrf.txt",list.files(i))]
    idf <- list.files(i, full.names = T)[grep(".idf.txt",list.files(i))]
    ls[["sdrf"]] <- read.delim(file = sdrf, row.names = 1)
    ls[["idf"]] <- read.delim(file = idf)
    for(k in AE$processedFiles) {
      ls[[k]] <- read.table(file = paste(i,k,sep = "/"),row.names = 1)
    }
    ls[["info"]] <- AE
  }
  setwd("../")
  if(import == T) return(ls)
}

# This function is ment for a bulk download without import into R. You can provide
# multiple accession IDs at once.
AEDownloadBulk <- function(accession, type = "processed", out = getwd()) {
  # create output dir
  dir.create(paste0(out,"/ArrayExpress"), showWarnings = F)
  setwd(paste0(out,"/ArrayExpress"))
  # run getAE() in loop for accession numbers
  for(i in accession) {
    dir.create(i, showWarnings = F)
    getAE(i, type = type, path = i, extract = T)
    zip <- list.files(i, full.names = T)[grep(".zip",list.files(i))]
    if(file.exists(zip[1])) file.remove(zip)
    }
  setwd("../")
}


## Now we use the function defined above to download & import the following accessions into R
accession = c("E-MTAB-9056","E-MTAB-9054")
for(i in accession) {
  x <- AEDownload(i)
  assign(i,x)
}
rm(x,i)

# If you only wish to download the files without loading them into R simply run: 
AEDownloadBulk(accession)

###################################


# Congrats! :) 
# Now you have successfully downloaded the processed data files for "E-MTAB-9056" and "E-MTAB-9054". 
# You should  have two R objects with those names in you R environment now.
# Search your environment for E-MTAB-9056" and "E-MTAB-9054" and look at their structures 
# (hint: str() and View())

# Q.4) What kind of objects are "E-MTAB-9056" and "E-MTAB-9054"? 
# How many objects are stored in them? (hint: str() and View())

# Q.5) Both E-MTAB objects should contain an IDF & SDRF file. What does IDF and SDRF stand for? What kind of information is stored in them?
# ( Hint: Check the vignette: browseVignettes("ArrayExpress") )


# Let's extract the count matrices from the E-MTAB objects
mtx.6PTU = `E-MTAB-9054`[["6PTU_CountMatrix.txt"]]
mtx.T3 = `E-MTAB-9056`[["T3_CountMatrix.txt"]]


# Q.6) Look at the CountMatrices. How many ENSEMBL Gene IDs and how many Samples (Sample IDs) are in there? 
# Do both experiments have the same number of samples? 


# Now let us have a look at the sample annotation files.
# Therefore we are extracting the sdrf file with the annotation info
anno.6PTU = `E-MTAB-9054`[["sdrf"]]
anno.T3 = `E-MTAB-9056`[["sdrf"]]

colnames(anno.6PTU)

# As you might see when you are looking at the, they contain a lot of columns.
# Many of those columns we don't need for downstream analysis. Hence only select the following:
coln <- c("Characteristics.growth.condition.","Comment.replicate.",
          "Comment.ENA_SAMPLE.","Characteristics.organism.","Assay.Name")
x = anno.6PTU[,coln]
x$Substance <- "6PTU"

coln2 <- c("Characteristics.growth.condition.","Comment.tank.spawning.group.",
           "Comment.ENA_SAMPLE.","Characteristics.organism.","Assay.Name")
y = anno.T3[,coln2]
y$Substance <- "T3"

# Replace column names
colnames(x)[c(1:2,5)] <- colnames(y)[c(1:2,5)] <- c("Condition","Tank","ID")

# Rename "normal" with "control" in the x$Condition column
x[x$Condition %in% "normal","Condition"] <- "control"

# Finally! Let us export our filtered annotation file to a csv
write.csv2(x, file = "ArrayExpress/E-MTAB-9054/6PTU_coldata.csv")
write.csv2(y, file = "ArrayExpress/E-MTAB-9056/T3_coldata.csv")


# Well done! :) Now you have everything you need for a differential gene expression analysis.
# The required input is a count matrix of genes and a sample annotation file (coldata file).
# From here on we can continue with EdgeR or DESeq2


## Session Information & save Rdata  ## ------------
sink(paste0(home,"/ArrayExpress/AE_SessionInfo.txt"))
print(date())
print(devtools::session_info())
sink()

save.image(paste0(home,"/ArrayExpress/AE.RData"))
# you can load this object via:
# load("path/to/file/biomaRt.RData"))

#######################################
###         END OF SCRIPT           ###
#######################################