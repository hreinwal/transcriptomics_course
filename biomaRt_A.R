###################
###   biomaRt   ###
###################

# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# README # ---------------------------------
# This is a training script to learn the basics of the biomaRt package.
# There is a great training Vignette available under:
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
##########

### Get the biomaRt package for R ### -----------------------
# Run these lines of code to install the required packages
# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("biomaRt")
#####################################

require(biomaRt)
vignette("biomaRt") # Check out the vignette
# Try the function browseVignettes() as well for biomaRt. What is the difference?

setwd("../") #go one dir up (out of the git repo)
home = getwd() #set your local working directory before you go to work


################################################
### Preparing a BioMart database and dataset ###
################################################

# Before we can use biomaRt for a search query to annotate our IDs, we first must connect to 
# an available ensembl server to retrieve the requiered and up to date information. 
# Luckily for us, the biomaRt package makes this task fairly easy.

?listMarts # check out the listMarts function. What is it doing?
listMarts(host = "https://www.ensembl.org")

?useMart # check out the useMart() function. How can I extract a dataset from it?
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")

?listDatasets # what is this function for?
datasets <- listDatasets(ensembl)

## Inspect the datasets object 
# 1) In which line can we find the zebrafish (danio rerio) dataset and how is it called?
#    Which version is the reference genome?
# Hint: try ?View() or ?searchDatasets()
# Solution -----------------
View(datasets)
searchDatasets(ensembl, pattern = "rerio")
# however much faster would be the search with grep():
datasets[grep("rerio",datasets$dataset),] #line 55
##############

# 2) Now that you know your dataset's name, store the d.rerio dataset 
# into an object called 'rerio'. (hint: ?useMart())
# Is it possible to do this without creating the ensembl object first?
# Solution ---------------------------
# various solutions
rerio = useDataset("drerio_gene_ensembl", mart=ensembl)
rerio = useDataset(datasets[55,1], mart=ensembl)
rerio = useDataset(datasets[grep("rerio",datasets$dataset),1], mart=ensembl)

# or simply in one step - most robust solution!
rerio = useMart("ENSEMBL_MART_ENSEMBL", dataset = 'drerio_gene_ensembl', host = "https://www.ensembl.org")
# or this. simple but works not always ... idk why to be honest ;) 
listEnsembl()
rerio = useEnsembl("genes", dataset = 'drerio_gene_ensembl')
##############

# Let's make a local copy of the rerio dataset mart-object for later usage 
# (in case server connection might fail later on)
save(rerio, file = paste0(home,"/drerio_mart.Robj"))

# Imagine you are a genius exotoxicologist working with Daphnia magna.
# 3) Can you find the Daphnia magna dataset in there? 
# 4) Does ensembl provide a D.magna reference genome? (Hint: listEnsemblGenomes() or search online)
# 5) How could we load the D.magna dataset using useMart()?
# Solution -----------------------------------
listMarts(host = "https://metazoa.ensembl.org")
magna = useMart("metazoa_mart", host = "https://metazoa.ensembl.org")
datasets2 <- listDatasets(magna)
datasets2[grep("magna",datasets2$dataset),] # line 51
magna = useDataset("dmagna_eg_gene", mart=magna)

# or in one run once you know what to look for:
magna = useMart("metazoa_mart", host = "https://metazoa.ensembl.org", 
                dataset = 'dmagna_eg_gene')
pulex = useMart("metazoa_mart", host = "https://metazoa.ensembl.org", 
                dataset = 'dpulex_eg_gene')

# alternatively for plant fans
ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
datasetsP <- listDatasets(ensembl_plants)
searchDatasets(ensembl_plants, pattern = "Arabidopsis")

# store a local copy of magna
#save(magna, file = "S:/data/biomaRt/dmagna_mart.Robj")
##############



###################################
###   Running a biomaRt query   ###
###################################

# Now that we know how to create a mart database object it's time to run a search query
# with the 'rerio' object we just created.

?getBM # check out this function. What is its purpose? What is the 'filters' parameter for?
# 6) Find a way to list all available filters for the rerio object (hint: ?listFilters())
# How many filters are there?
# Solution -----------------------
filters = listFilters(rerio) 
nrow(filters) #410
View(filters)
#############

# 7) What is the 'attributes' parameter for? How many 'attributes' can we apply on the 'rerio' mart object?
# Solution -----------------------
attrib = listAttributes(rerio)
nrow(attrib) #2892
View(attrib)
##############

## Side note:
# If you are searching for a particular attribute or filter you can use:
# ?searchAttributes()
# ?searchFilters()



#################################
### Annotate ensembl gene IDs ###
#################################

# In this next section we will practice the search query via biomaRt to learn how
# to effectively annotate ENSEMBL Gene IDs.

# This is a vector list of IDs we wish to annotate for now:
ID <- c("ENSDARG00000079305","ENSDARG00000045142","ENSDARG00000059139",
        "ENSDARG00000095767","ENSDARG00000053481","ENSDARG00000092660",
        "ENSDARG00000001975","ENSDARG00000099351","ENSDARG00000053136",
        "ENSDARG00000054588","ENSDARG00000117519","ENSDARG00000082287",
        "ENSDARG00000081280","ENSDARG00000087732","ENSDARG00000057055")

# ENSEMBL gene id is the type of ID we use here
idType <- "ensembl_gene_id"

# The information we would like to retrieve 
attr <- c("ensembl_gene_id","ensembl_peptide_id","external_gene_name",
          "description","gene_biotype","zfin_id_id",
          "entrezgene_accession","entrezgene_id")
# you can also retrieve: "reactome", "go_id", "kegg_enzyme" attributes 

# Now that we have all the required parameters fixed for the query,
# let's annotate the IDs
anno <- getBM(attributes = attr, mart = rerio, uniqueRows = T,
              filters = idType, 
              values = ID)
View(anno)

# 8) Inspect 'anno'. (View()) Is the length of 'anno' (nrow()) the same as 
#   the length of ID? Which columns contain redundant information?

# 9) Which attribute produces duplicated results and why? How could we avoid that?
# Solution ---------------------------
# ensembl_peptide_id attribute creates duplicates due to splicing variants and poor annotation
attr1 <- c("ensembl_gene_id","external_gene_name","description","gene_biotype",
          "zfin_id_id","zfin_id_symbol")
anno1 <- getBM(attributes = attr1, mart = rerio, uniqueRows = T,
                       filters = idType, 
                       values = ID)
View(anno1)
##############

# 10) Are all of the annotated ENSEMBL IDs 'idType' protein coding sequences?
#     How could you retrieve only "protein_coding" IDs? Is the biotype annotation accurate?
#     Hint: try to use the 'filters' parameter in the getBM() function. ?getBM()
# Solution -----------------
# Annotate only the protein coding genes
anno.Prot <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"), 
                   mart = rerio, uniqueRows = T,
                   filters = c(idType,"biotype"), values = list(ID,"protein_coding"))
View(anno.Prot)
##############



######################################
###  Annotating our DESeq2 results ###
######################################

# Now let's load DESeq2 result tables we created the previous days and annotated these tables.
# We will use the log2-fold shrinked result tables for  that.
filePath <- "DESeq2_Pairwise/Results/" #provide here the path to your DESeq2 results

tmp = list.files(filePath, pattern = "_reslfs_", full.names = T)
T3_HE <- read.csv2(tmp[1], row.names = 1)
T3_LE <- read.csv2(tmp[2], row.names = 1)
# Have a look at the data via:
head(T3_HE)
# Note: The row names of these data tables are the ENSEMBL Gene IDs

#  Now try on your own to create an annotation object called 'anno2' for all ENSEMBL gene IDs with the following attributes:
# "ensembl_gene_id","external_gene_name","description","gene_biotype","entrezgene_accession"
#  You can select the ID you want to annotate i.e. via rownames(T3_HE)
# 
# Solution ---------------------
id     <- rownames(T3_HE)
idType <- "ensembl_gene_id"
attr2  <- c("ensembl_gene_id","external_gene_name","description","gene_biotype")
# In case you need to recreacte the mart object 'rerio' again: 
# rerio <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'drerio_gene_ensembl', host = "https://www.ensembl.org")
anno2  <- getBM(attributes = attr2, mart = rerio, uniqueRows = T,
              filters = idType, values = id)
##############

# Using ?merge() we can simply bring together these two data tables. The easiest way is merge
# them by rownames. Hence we first name the rownames of 'anno2' after the anno2$ensembl_gene_id column
str(row.names(anno2))
row.names(anno2) <- anno2$ensembl_gene_id
str(row.names(anno2))
# Quick check if all IDs are there
stopifnot(length(row.names(anno2)) == length(row.names(T3_HE)) &
            length(row.names(anno2)) == length(row.names(T3_LE)))

# Now we merge the files:
T3_HE.an <- merge(T3_HE, anno2, by  = 0)[,-1]
T3_LE.an <- merge(T3_LE, anno2, by  = 0)[,-1]
# And add row names again
row.names(T3_HE.an) <- T3_HE.an$ensembl_gene_id
row.names(T3_LE.an) <- T3_LE.an$ensembl_gene_id

# Done! :) Great job! You successfully annotated your DEseq2 results with meaningful biological data.
# This enables us now to further explore the biological functions of the genes observed.

# Let's have a look of the different gene biotypes represented in our dataset. 
x = table(T3_HE.an$gene_biotype)
x
# 11) How many gene biotypes are in your dataset? Which is the largest set with how many counts?
# Solution --------------
length(x) #24
which.max(x) #protein_coding 
max(x) #22740
##############



#########################
###   Data plotting   ###
#########################

# A very easy, quick and dirty way to visualize the biotype data:
barplot(x, main = "Gene biotypes", ylab = "Counts")
barplot(sort(x, decreasing = T), main = "Gene biotypes", ylab = "Counts")

# Much nicer with ggplot2:
require(ggplot2)
require(ggpubr)
ggbar <- function(x, subtitle="", xlab =""){
  # create df for plotting
  df <- as.data.frame(x)
  df$Var1 <- factor(df$Var1, levels = names(sort(x, decreasing = T)))
  S <- sum(df$Freq)
  df$perc <- round((df$Freq/S)*100,2)
  colnames(df)[1] <- "Biotype"
  # plot
  ggplot(df, aes(x = Biotype, y = Freq, fill = Biotype)) + theme_light() +
    theme(panel.grid.minor = element_blank(), legend.position = "none") +
    labs(subtitle = subtitle, x = xlab, y = 'Counts') + rotate_x_text(angle = 25) +
    theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13)) +
    geom_bar(stat="identity", alpha=.8, width = .8) +
    geom_text(aes(x=Biotype, y=Freq+(max(Freq)*.03), label= paste0(perc," %")), angle = "0")
}
ggbar(x)
# only the biotypes with a count greater 5
ggbar(x[which(x > 5)])



## Session Information & save Rdata  ## ------------
sink(paste0(home,"/biomaRt_SessionInfo.txt"))
print(date())
print(devtools::session_info())
sink()

save.image(paste0(home,"/biomaRt.RData"))
# you can load this object via:
# load("path/to/file/biomaRt.RData"))
#######################################

#############################################
###   Yay well done! Done for today! :)   ###
#############################################


### Additional questions in case the stuff above was to easy ;) ###

# Supl.Q1: Can you find a way to filter 'T3_HE.an' and 'T3_LE.an' only for:
# protein_coding, processed_transcript, processed_pseudogene and  polymorphic_pseudogene?
# Hint: look into ?subset() or the %in% filter option.
# Solution: ---------------------------
df <- T3_HE.an
val <- c("protein_coding","processed_transcript","processed_pseudogene","polymorphic_pseudogene")
res <- df[df$gene_biotype %in% val,]
ggbar(table(res$gene_biotype))
###############

# Supl. Q2: Match ortholog genes between species using ?getLDS()
?getLDS



####  END OF SCRIPT  ####