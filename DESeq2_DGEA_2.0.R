#####################################################################
### DESeq2 RNA Seq DGE analysis for multiple factors / treatments ###
#####################################################################

# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# README # ----------------------------------------------------------------------
# This script will perform data norm and differential gene expression analysis on 
# RNASeq count data with the DESEq2 package. It will perform some basic quality measure
# plottings (variance distr., PCA, dissimilarity matrix, pvalue distr.,...)

# Requiered Input: 
# - CountMatrix 
# - coldata.csv

# As it is easier to format the coldata file than the count matrix,
# the samples / conditions listed in your coldata file will be tested.
# If you wish to exclude certain parameters from the analyis just remove them 
# from the coldata file prior running the sript.

# Main Analysis Steps are: 
# 1) RLE normalization (DESeq2) across all samples element of CountMatrix
# 2) Check Log2FC distribution (with respect to control) and determine LFcutOff 
#    (Biological effect size cutoff) upper 90% quantile of abs(LFC values))
# 3) apeglm shrinking on LFC (log2FC shrinkage)
# 4) Multiple t-testing with Benjamin-Hochberg correction (padj < 0.05) 
#    and independent hypothesis weighing (IHW) to identify DEGs
#    for LFC and apeglm(LFC) values for H0: LFC = 0; apeglm(LFC) = 0
##########

### LOAD PACKAGES ### -------------------------------------------------
require(RColorBrewer)
require(ggplot2)
require(ggpubr)
require(pheatmap)
#require(hexbin)
require(DESeq2) #installed
require(IHW)    #installed
require(apeglm) #installed
require(qvalue)
#require(locfdr)
#####################

# Navigate to the directory containing your countmatrix and coldata file.
setwd("../data_T3/")

# Makes sure that R will find your files by checking via:
pathToData = list.files(pattern = "coldata[.]csv|[Cc]ount[Mm]atrix[.]", recursive = T, full.names = T)
# Look at the pathToData object you just created.
pathToData
# R should display you two file paths. One for coldata and the other for the count matrix.
# If not you are in the wronk working directory. Please check again ;)
# If R displays more than 2 files you most likely generated more than one count matrix yesterday.
# Please remove duplcated ones before you continue.
stopifnot(length(pathToData) == 2) # <--- can you figure out what this line is for?

# Q: Which input files are required to run DESeq2 and what info do they contain?


### Import data ### ---------------------------------------------
countMtx = read.table(grep(pathToData, pattern = "[Cc]ount[Mm]atrix", value = T), row.names = 1)
coldata = read.csv2(grep(pathToData, pattern = "coldata[.]", value = T), row.names = 1)
# Q: You should have an object called countMtx & coldata in your R environment now. Check their structure
# using the dim() and str() function. If you are not sure what they do - remember you can use the ? to get help.
# e.g. type: ?str()


## format input files -------------------------------------------
# Q: check the structure of your coldata object. Are there any factors assigned yet?
str(coldata)

# set factor levels
coldata$Tank      <- factor(coldata$Tank)
coldata$Condition <- factor(coldata$Condition)
coldata$Substance <- factor(coldata$Substance)
# Q: now check structure of coldata again. What changed? 

#Extract the name(s) of the tested substance(s)
substance <- levels(coldata$Substance)

# relevel factors in Condition in right order:
condition = levels(coldata$Condition)

# Q: What is the order of the levels in condition ? Are you ok with the order of factors?
# to check for that run: levels(coldata$Condition)
# To rearrange the factor level order run the next lines....
tmp = as.character(coldata$Condition[grep("[Cc]ontrol",coldata$Condition)][1]) #extract control factor level
if(all(grepl("[Ee]xposure",condition[2:4]))==T & length(condition)==4) {# Control, LE, ME, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,4,2)])
} else if (all(grepl("[Ee]xposure",condition[2:3]))==T & length(condition)==3) { # Control, LE, HE
  coldata$Condition <- factor(coldata$Condition, levels = condition[c(1,3,2)])
} else { # just specifying the reference level --> Control
  coldata$Condition <- relevel(coldata$Condition, ref = tmp)
}
message(paste0("Reference level:\t",levels(coldata$Condition)[1]),"\nSorting ",length(condition)-1," Treatments:\t",paste(levels(coldata$Condition)[-1], collapse = " "))
condition = levels(coldata$Condition)
# Q: What is the order in condition now? > levels(coldata$Condition)


# Sort & check correct data import:
countMtx <- countMtx[,rownames(coldata)] # reduce the to analyse samples in mtx by the samples listed in coldata
stopifnot(all(rownames(coldata) == colnames(countMtx))) # Must be TRUE!

# Proper data formating is a crucial step in R data analysis. This usually takes up most of the time.
# Once you have everything properly imported and setup downstream steps are much easier to perform.

rm(tmp, pathToData) # <--- can you figure out what this command does?
###################



##############################################################################################
###### All imported! Yay! :) Let's START DESeq2 differential gene expression analysis! #######
##############################################################################################


### Setup output environment ###
dir.create("DESeq2_Pairwise", showWarnings = F)
setwd("DESeq2_Pairwise")
home=getwd()

p = .05 # padj cut off for any downstream analysis

##################
###   DESeq2   ###
##################

browseVignettes("DESeq2") # this will take you to the official DESeq2 vignette which explains this
# package and its powerful functions and features in great depth. If you ever come back using
# DESeq2 I highly recomend to read over it. (There are also a lot of very good code examples and
# coding recommendations). Helped me a lot for my data analysis.

## Import data into a single DESeq2  object ## -------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countMtx, # <--- our count matrix
                              colData = coldata,    # <--- our coldata file
                              design = ~ Tank + Condition) # Multifactor Level design to correct for batch effect

# You successfully created now a DESeq2 object in which the combined information of 
# coldata and count matrix as well as our test design is stored. By default DEseq2 takes the
# first level of the factor (in our case condition) which we wish to compare for as reference level

# To access i.e. the count matrix stored in dds you simply use counts()
counts(dds) # prints the same object as countMtx! :) <--- what would be a better way to display this?


## Count Matrix filtering ## ---------------------------------------
# Before we start with statistical test for DGE analysis we should remove low read counts
# from our count matrix. But how do we determine a "good" threshold? 
# The truth is ... there is NO standard procedure for that. A very common method is to 
# take remove every gene which has a rowSum smaller than a certain value. (very common 10 - 20)

# However, this cutoff by an arbitrary value is biased simply by the total number of samples
# you have in your count matrix. (More samples -> easier to reach the threshold)

# Therefore it is a good idea to look at the quantile distribution to get an idea 
# of a reasonable cutoff size.
# computing rowSum for all rows: 
tmp = apply(countMtx, 1, sum) # computes the sum of rows
tmp1 = quantile(tmp, probs = seq(0, 1, 0.05)) # Here we look at the 5% quantile invervals from 0 to 100%
tmp1
# Q: Which quantile contains more than 0 counts?


# We can simply plot the results with:
par(mfrow=c(1,2)) # creates a frame to show plots in with 1 row & 2 columns
df <- data.frame(x=seq(0,1,.05)*100, y=tmp1)
plot(df, xlab = "Quantile-%", ylab = "Total counts", main = "RowSums")

# This doesn't show us to much does it? The largest values are just to large to make
# any good visual comparison. The trick is to transform the data to bring the values 
# closer together. We will run a simple log10 transformation here. 
# Keep in mind the log10(0) is not defined! => we add the value of 1 to each count in the matrix
tmp2 = quantile(log10(tmp+1), probs = seq(0, 1, .05))
df <- data.frame(x=seq(0,1,.05)*100, y=tmp2)
plot(df, xlab = "Quantile-%", ylab = "log10(N+1)", main = "transformed rowSums")
par(mfrow=c(1,1))
# This looks much nicer doesn't it? :)

# Removing the lower 20 - 30 % of your counts is usually a reasonable approach.
# I made the experience that the number of samples in your count matrix (ncol)
# is usually a good threshold for that which scales with your sample size. 

# For now let us use a minimum rowSum threshold of NumberOfSamples + 2
# Q: What is our rowSum threshold then? 

# How the filtering affects the number of genes analyzed with DESeq2 
# can be visualized via:
cutOffPlot <- function(countMtx, cut = ncol(countMtx)) {
  n = ncol(countMtx)
  if(n < 6) {stop("Number of count matrix columns < 6")}
  X = seq(n-6,n+50,1) 
  X[X == 0] <- 1
  X = sort(union(X, c(1:5)))
  
  rNbr = c() # empty vector to store Nbr of genes per cutoff in
  for(x in X) {
    tmp  = nrow(countMtx[which(rowSums(countMtx) >= x), ])
    rNbr = append(rNbr, tmp)
  }
  df <- data.frame(x = X, y = rNbr)
  plot(df, xlab = "Row sum cut off", ylab = "Nbr of genes",
       main = "Genes in count matrix")
  abline(v=cut, lwd=1.5, lty=2, col = "firebrick")
  #abline(h=df[cut,2], lwd=1.5, lty=2, col = "blue")
  text(x = cut, y = df[cut,2], pos = 4, offset = 1.5,
       labels = paste("Cutoff:",cut,"=",df[cut,2],"genes"))
}
par(mfrow=c(3,1))
cutOffPlot(countMtx)
cutOffPlot(countMtx, 11)
cutOffPlot(countMtx, 20)
par(mfrow=c(1,1))

# Filter count matrix - removing genes with none or low counts
# selecting gene names / rownames from count matrix for which rowSum is > ncol +2
keep <- rowSums(counts(dds)) >= ncol(counts(dds))+2
dds  <- dds[keep,] #subsetting dds with selection of tmp
# Q: how many genes are now left? (hint: nrow() counts())

## This is actually a good time now to check the distribution of your read counts
# Let's compute the row mean values for the filtered and unfiltered count mtx
tmp = apply(countMtx, 1, mean)
tmp1 = apply(countMtx[keep,], 1, mean)

# plotting
par(mfrow = c(1,3))
hist(log2(tmp+1), breaks = 30, main = "countMtx", 
     ylim = c(0,2500), xlim = c(0,20),
     xlab = "log2(Row mean count +1)")
hist(log2(tmp1 +1), breaks = 30, main = "filtered countMtx",
     ylim = c(0,2500), xlim = c(0,20),
     xlab ="log2(Row mean count +1)") #histogram of filtered counts
cutOffPlot(countMtx, 11)


## Running DESeq2 & first data QC ## ----------------------------------------------------------------
# Now we can finally run DESeq2 :) Yay!
# use ?DESeq to inspect the function.
# Q: Which main steps are performed by this function? 
# Q: Which model ("fitType") is the default?

message(paste0("\n Starting DESeq2 Analysis - Pairwise Wald's t-test and IHW \n padj < ",p,"\n Tested Substance: \t",substance,"\n"))
dds <- DESeq(dds, test = "Wald")                   # Pairwise comparison


# If you wish to run the ANOVA like LRT method you can run it like that:
#dds <- DESeq(dds, test = "LRT", reduced = ~ Tank) # ANOVA-like approach

# That was quick wasn't it? :) Congrats! You just finished DGE analysis with DESeq2.
# Well ... partially ... we still need to extract the results and plot them ;) 

# You can quickly check if your normalization has worked by plotting: 
par(mfrow=c(2,2))
barplot(colSums(counts(dds)),
        ylab= "Total gene counts",
        main= "Raw counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
barplot(colSums(counts(dds, normalized=T)),  #plot mean normalized counts
        ylab= "Total gene counts",
        main= "DESeq2 norm. counts",
        las= 3, #rotating sample labels 90
        ylim= c(0,1.2*max(colSums(counts(dds)))))
boxplot(log10(counts(dds)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "Raw counts")
boxplot(log10(counts(dds, normalized=T)+1),
        ylab= "Log10(Gene counts + 1)",
        las = 3, #rotating sample labels 90
        main= "DESeq2 norm. counts")

# Q: Do you think the moralization has worked?
##############################################


### Extracting DESeq2 results - normalized count matrix ### -----------------------------------------------------
# In the first step you might want to extract the normalized count matrix. 
# This is fairly easy using the counts() function which you already know.

# Normalized raw count matrix
normMtx <- round(counts(dds, normalized = T),3) 

# Variance stabilizing transformation on normalized count matrix
vst.bl  <- assay(vst(dds, blind = T)) # blind=T; use for QC. Assesses data unbiased by Condition or Tank
vst     <- assay(vst(dds, blind = F)) 

# log2 transformed normalized count matrix
ntd <- assay(normTransform(dds)) # Extract (n+1)log2 transformed mean read counts


# From all the matrix objects we created above,
# vst.bl & vst are the most suitable for PCA, t-SNE, dissimilarity Mtx, heatmap etc ....

# How these different data transformation methods affect your data you can 
# check with a simple boxplot:
par(mfrow=c(1,4))
boxplot(normMtx, notch = TRUE,
        las = 3, #rotating sample labels 90
        main = "Normalized read counts",
        ylab = "norm. read counts")

boxplot(ntd, notch = TRUE,
        las = 3, #rotating sample labels 90
        main = "log2 Transformation",
        ylab = "log2 (norm. read counts + 1)")

boxplot(vst.bl, notch = TRUE,
        las = 3, #rotating sample labels 90
        main = "vst.bl Transformation",
        ylab = "vst (norm. read counts)")

boxplot(vst, notch = TRUE,
        las = 3, #rotating sample labels 90
        main = "vst Transformation",
        ylab = "vst (norm. read counts)")
###########################################################


### Sample correlation ### --------------------------------------------------------------
# Biological sample correlation plots can help you to identify if i.e. biological replicates show
# a good level of correlation. This means we assess if our samples "behave" in a similar fashion.
# If two samples show similar expression profiles we would expect to see a high level of correlation
# (at least for the biological replicates)

# here we are building a custom plot function
corplot.1 <- function(df1) { ggplot(df1) +
  geom_point(aes(x, y, col = d), size = .8, alpha = .4) +
  labs(x=paste0(ID[p1],transf),y=paste0(ID[p2],transf),
       title = (paste0(substance,": ",i)),
       subtitle = (paste("[",coldata[ID[p1],"Tank"],"vs",coldata[ID[p2],"Tank"],"]"))) +
  annotate("text",
           label = paste0("Pearson = ",round(cor(df1$x,df1$y),2)),
           x = (max(df1$x)), 
           y = (min(df1$y)),
           hjust=1, vjust=0) + 
  annotate("text",
           label = paste0("R2 = ",round((cor(df1$x,df1$y))^2,2)),
           x = (min(df1$x)), 
           y = (max(df1$y)),
           hjust=0, vjust=1) + 
  scale_color_identity() +
  coord_equal(ratio=1) + 
  theme_bw() }

# Log10 transformed mean counts
df      <- ntd #INPUT 
transf  <- " - [log2(counts+1)]"
gg.list <- list() # empty list to store outputs in
for (i in condition){
  ID <- rownames(coldata)[which(coldata$Condition %in% i)]
  p1 <- 1 #ID vector position
  p2 <- 2 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 1 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 2 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
}
# Export to png
png(filename = paste0(substance,"_Correlation_log10.png"),
    width = 7.33333*3, height = 7.6*length(condition), units = "cm", bg = "white",
    pointsize = 1, res = 450)
print(ggpubr::ggarrange(plotlist = gg.list, ncol = 3, nrow = (length(ls(gg.list))/3)))
dev.off()


# vst transformed mean counts
df      <- vst #INPUT
transf  <- " - [vst(counts)]"
gg.list <- list() # empty list to store outputs in
for (i in condition){
  ID <- rownames(coldata)[which(coldata$Condition %in% i)]
  p1 <- 1 #ID vector position
  p2 <- 2 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 1 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
  
  p1 <- 2 #ID vector position
  p2 <- 3 #ID vector position
  df1 <- data.frame(x = df[,ID[p1]],
                    y = df[,ID[p2]],
                    d = densCols(df[,ID[p1]], df[,ID[p2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  )
  gg.list[[paste0(i,".",ID[p1],"vs",ID[p2])]] <- corplot.1(df1)
}
# Export to png
png(filename = paste0(substance,"_Correlation_vst.png"),
    width = 7.33333*3, height = 7.6*length(condition), units = "cm", bg = "white",
    pointsize = 1, res = 450)
print(ggpubr::ggarrange(plotlist = gg.list, ncol = 3, nrow = (length(ls(gg.list))/3)))
dev.off()

rm(df,ggobj,transf,gg.list,i,p1,p2,df1)

# Q: The plots were exported. Have a look at them. Which samples are most similiar?
##########################


### PCA ### ---------------------------------------------------------
require(pcaMethods) # we will use this package for the PCA
# here is a small function you can use to perform and quickly plot a PCA.
# the only things you need to provide is the value matrix and a sample annotation file (coldata)
myPCA <- function(mtx, coldat, pcaM ="svd", top = 2000, title = "") {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting PCA for Var.Top:\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # running pca() - this line is all it takes for a PCA
  pc <- pca(X, method = pcaM, center = T, nPcs=2)
  pcaDf <- merge(coldat, scores(pc), by=0)
  
  # plotting color
  col.Cond = colorRampPalette(brewer.pal(n=8, name="YlOrRd"))(length(levels(coldat$Condition)))
  col.Cond[1] <- "#56B4E9" #Defines color for control
  names(col.Cond) <- levels(coldat$Condition)
  
  # plotting the PCA results
  ggplot(pcaDf, aes(PC1, PC2, colour = Condition, shape = Tank)) +
    geom_point(size = 3, alpha = .65) + scale_colour_manual(values = col.Cond) +
    ggtitle(paste0(title," - expl.Var[%]: ",round(pc@R2cum[2]*100,1)," - ",pcaM," - Top:",pc@nVar)) +
    xlab(paste0("PC1: ",round((pc@R2)[1]*100,1),"% variance")) +
    ylab(paste0("PC2: ",round((pc@R2)[2]*100,1),"% variance")) + #stat_ellipse() +
    theme_bw() + theme(aspect.ratio = 1)
}

myPCA(vst, coldata, title = "vst")
myPCA(vst.bl, coldata, title = "vst.bl")
# Q: Compare the two PCA plots. Which one is "better"? Why are they different?
# Q: The "top" parameter refers to the topmost varying genes (largest variance among groups).
#    What happens when you change this parameter? (i.e. try 100 500 10000 20000)
###########


### Dissimilarity Matrix ### ----------------------------------------
myDismtx <- function(mtx, coldat, method = "euclidean", 
                     top = 2000, title = "", ...) {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting Sample Dist for Var.Top:\t\t",top))
  
  # annotation colors
  col.Cond = colorRampPalette(brewer.pal(n=8, name="YlOrRd"))(length(levels(coldat$Condition)))
  col.Cond[1] <- "#56B4E9" #Defines color for control
  col.Tank = colorRampPalette(c("gray95","gray50"))(length(levels(coldat$Tank)))
  
  ann_colors = list(Tank = setNames(col.Tank, levels(coldat$Tank)),
                    Condition = setNames(col.Cond, levels(coldat$Condition)))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # transpose input, calculate sample distance and create a distance matrix
  sampleDist <- dist(X, method = method) #must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  distMtx <- as.matrix(sampleDist)
  rownames(distMtx) <- paste(coldat$Condition, coldat$Tank, sep="-")
  
  # create annotation object for heatmap
  ann_col <- subset(coldat, select = "Condition")
  ann_row <- subset(coldat, select = "Tank")
  rownames(ann_row) <- paste(coldat$Condition, coldat$Tank, sep="-")
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(15)
  
  # plot
  pheatmap::pheatmap(
    distMtx, angle_col = "90", drop_levels = T,
    number_format = if(method == "manhattan"){"%1.0f"}else{"%.1f"},
    main = paste0(title," - SampleDist [",method,"] - Var.Top: ",top),
    clustering_distance_rows = sampleDist,
    clustering_distance_cols = sampleDist,
    annotation_col = ann_col, # Condition!
    annotation_row = ann_row, # Tank!
    annotation_colors = ann_colors,
    col = colors,
    ...)
}

myDismtx(vst, coldata, title = "vst transformed norm counts")
myDismtx(vst.bl, coldata, title = "vst.bl transformed norm counts")

# Play around with the parameters a little bit and see how the clustering changes.
myDismtx(vst.bl, coldata, title = "vst.bl transformed norm counts", top = 10000,
         clustering_method = "average")
# Q: What is the difference between average, complete and single linkage in hclustering?

############################


### Extracting DESeq2 results - differential expressed gene tables ### ------------------------------------
# With resultsNames you can inspect the names of the different pairwise comparisons DESeq2 performed
resultsNames(dds)
# As you might remember we are only interested in the Condition comparison.
tmp = resultsNames(dds)[grepl("^[Cc]ondition_",resultsNames(dds))]
tmp # these are the objects we will now retrieve from the large DESeq2 object
# for this purpose DESeq2 provides the results() function. With it we can retrieve the results
# Check: ?results()
# Let's save the results in a list object 
res.ls <- list() #create an empty list to store objects in
for(i in tmp) {
  x = gsub("[Cc]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, non-shrunk Log2FC]"))
  res.ls[[gsub("_.+","",x)]] <- results(dds, name = i, filterFun = ihw, alpha = p)
  message("Done!\n")
}

# You can simply access the object in the list via the $ extension:
res.ls$LowExposure
# with summary() you can get some quick stats
summary(res.ls$LowExposure)
# Q: How many genes are significantly upregulated? How many downregulated? How many outliers?
# Check that for High and Low exposure condition.

# To check the distribution of pvalues you can run:
par(mfrow=c(1,2))
hist(res.ls$LowExposure[["pvalue"]], main= "pval distr. - Low treatment",
     col = "gray50", border = "gray50", ylab = "Frequency", xlab = "pvalue",
     breaks = 100)
# and for high exposure treatment: 
hist(res.ls$HighExposure[["pvalue"]], main= "pval distr. - High treatment",
     col = "gray50", border = "gray50", ylab = "Frequency", xlab = "pvalue",
     breaks = 100)

# Or you can run it in a loop: 
par(mfrow=c(1,length(res.ls)))
for(i in names(res.ls)) {
  hist(res.ls[[i]][["pvalue"]], main = paste("pval distr.",i),
       col = "gray50", border = "gray50", ylab = "Frequency", xlab = "pvalue",
       breaks = 100)
}
par(mfrow=c(1,1))


# DESeq2 also provides the ?lfcShrink() function wich implements the so called
# apeglm shrinkage method. This method was published in 2018 by Zhu et al. and
# allows a powerful transformation of the log2FC values in our deseq2 result table.
# As stated in the name, this transformation performs a shrinking of the Log2FC values.
# How much the Log2FC values are shrinked depends on how statistically reliable the value is.
# Hereby apeglm considers the overall baseMean of a particualr log2FC and its SD.
# Large SD values will be more penalized compared to log2FC values with low SD.
# Same goes for genes with a low baseMean count. Hence low abundant genes with a large
# SD will be heavily penalized. This means their log2FC value will be heavily shrinked!
# Note: We observed that (unless you apply a particular)
resLfs.ls <- list() # save log2FC shrunk results in different list
for(i in tmp) {
  x = gsub("[cC]ondition_","", i)
  message(paste0("Extracting DESeq2 results for: ",x," [IHW, apeglm shrunk Log2FC]"))
  resLfs.ls[[gsub("_.+","",x)]] <- lfcShrink(dds, coef = i, type = "apeglm", 
                                             res = res.ls[[gsub("_.+","",x)]])
  message("Done!\n")
}
rm(x,i)
# The shrunk log2FC values are in particular useful for data visualization such as MA- & Vulcano plots
# Additionally they are better for downstream Gene Set enrichment analysis as the data is less "noisy".
######################################################################


### Annotating ENSEMBL Gene IDs with BiomaRt ### ----------------------------------------------------
require(biomaRt)
# biomaRt is a package designed for effectively searching ENSEMBL database for gene annotation.
# Please check out the biomaRt.R script to learn how to use it. If you are planing to work
# with large set genomic data file I highly recommend you to get familiar with this package
# as it will make your life a whole lot easier. browse the vignette here: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

# For now it is ok if you just run the next lines of code which will automatically annotate
# your DESeq2 result tables.

rerio  <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'drerio_gene_ensembl', host = "https://www.ensembl.org")
id     <- rownames(res.ls[[1]])
idType <- "ensembl_gene_id"
attr  <- c("ensembl_gene_id","external_gene_name","description","gene_biotype"#,"entrezgene_id"
           )
GeneAnno <- getBM(attributes = attr, mart = rerio, uniqueRows = T,
                filters = idType, values = id, useCache = F )
row.names(GeneAnno) <- GeneAnno$ensembl_gene_id

# Now merge GeneAnno with DESeq2 result tables
AnnoFun <- function(x){
  tmp = merge(as.data.frame(x), GeneAnno, by  = 0)[,-1]
  row.names(tmp) = tmp$ensembl_gene_id
  tmp[,-which(colnames(tmp) %in% "ensembl_gene_id")]
}
res.ls <- lapply(res.ls, FUN = AnnoFun)
resLfs.ls <- lapply(resLfs.ls, FUN = AnnoFun)

################################################


### Exporting result to csv ### -------------------------------------------------------
# Before we finally get to plot our results let us export the DESeq2 result tables to a csv file (Xcel readable)
dir.create("Results", showWarnings = F)
for(i in names(res.ls)){
  n = gsub("[Ee]xposure","",i)
  message(paste("Saving DESeq2 result table",i,"[IHW, non-shrunk Log2FC] to csv ..."))
  write.csv2(res.ls[[i]],
             file = paste0("Results/",substance,"_res_",n,".csv"))
  message(paste("Saving DESeq2 result table",i,"[IHW, apeglm shrunk Log2FC] to csv ..."))
  write.csv2(resLfs.ls[[i]],
             file = paste0("Results/",substance,"_reslfs_",n,".csv"))
  message("Done!\n")
}
rm(n,i)
###############################


### Vulcano and MA plots ### --------------------------------
# Here are two ready to use functions to MA & vulcano plot your DESeq2 result tables
MAfun <- function(res, LFcut, topgenes = 10, title = "", Symbol = F, shrink = F){
  if(missing(LFcut)){
    LFcut = quantile(abs(res$log2FoldChange), na.rm=T, .9) # 90% quantile of the abs(log2FC)
  }
  if(LFcut > 1){LFcut = 1} # in case the LFcut value is > 1, set LFcut to 1
  
  ggpubr::ggmaplot(
    res, 
    fdr = p, fc = 2^(LFcut), size = 1.5,
    top = topgenes,
    select.top.method = 'fc', # fc or padj
    palette = c("#B31B21", "#1465AC", "darkgray"),
    main = paste(substance,title,"[ LFcut:",round(LFcut,2),"]"),
    legend = "bottom",
    genenames = if(Symbol == F){NULL}else{as.vector(res$external_gene_name)},
    ylab = if(shrink == F){bquote(~Log[2]~ "fold change")}else{bquote("apeglm ("~Log[2]~"fc)")},
    xlab = bquote(~Log[2]~ "mean expression"),
    font.label = c("bold", 11), label.rectangle = F,
    font.legend = "bold",
    font.main = "bold",
    ggtheme = ggplot2::theme_light()
  )
}
VulcFun <- function(res, LFcut, topgenes = 10, title = "", Symbol = F, shrink = F,
                    pcut = .05){
  
  if(missing(LFcut)){
    LFcut = quantile(abs(res$log2FoldChange), na.rm=T, .9) # 90% quantile of the abs(log2FC)
  }
  if(LFcut > 1){LFcut = 1}
  
  if(Symbol == T){
    select <- res[order(res$padj),"external_gene_name"]
    } else {
      select <- rownames(res[order(res$padj),])
      }
  
  #if(shrink == F){bquote(~Log[2]~ "fold change")}else{bquote("apeglm ("~Log[2]~"fc)")}
  EnhancedVolcano::EnhancedVolcano(
    res, x = "log2FoldChange", y = "padj",
    title = paste(substance,title,"[ LFcut:",round(LFcut,2),"]"),
    subtitle = NULL,
    lab = if(Symbol == F){rownames(res)}else{res$external_gene_name},
    selectLab = if(topgenes == 0){select[NULL]}else{select[1:topgenes]}, # select topgenes based on padj value
    #legendLabels = c('NS','Log2 FC','padj','padj & Log2 FC'),
    xlab =  if(shrink == F){bquote(~Log[2]~ "fold change")}else{bquote("apeglm ("~Log[2]~"fc)")},
    ylab = bquote(~-Log[10]~italic(padj)),
    FCcutoff = LFcut,  #default 1
    pCutoff = pcut,       #default 10e-6
    labSize = 3.0,
    pointSize = 1.5, #default 0.8
    col = c("grey30", "grey30", "royalblue", "red2"),
    shape = c(1, 0, 17, 19),   #default 19, http://sape.inf.usi.ch/quick-reference/ggplot2/shape for details
    colAlpha = 0.4, #default 0.5
    hline = c(0.01, 0.001),
    hlineCol = c('grey40','grey55'),
    hlineType = 'dotted',
    hlineWidth = 0.6,
    gridlines.major = T,
    gridlines.minor = T,
    drawConnectors = TRUE,
  ) + theme_light() + theme(legend.position = "none")
}

# Now we can plot all results from our res.ls with a loop
gg <- list() #empty list to store our objects to plot in with
for(i in names(res.ls)) {
  res = res.ls[[i]] # picking the i result table
  gg[[paste0("MA.",i)]] <- MAfun(res, Symbol = T, title = i, topgenes = 6)
  gg[[paste0("VO.",i)]] <- VulcFun(res, Symbol = T, title = i, topgenes = 6)
}
# Now you can plot each plot of the gg list object individually. i.e:
names(gg)
gg[[1]]
gg[[2]]

# Or you can plot all plots in one panel:
ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = length(res.ls))

# It is the best to export this nice graph into a png.
# Remember pdf is better but in this case the multiple layers in the pdf might create 
# some problems for so many data points!
png(paste0(substance,"_MA_Volcano.png"),units = "cm", bg = "white", res = 500,
    width = 7.33333*3, height = 7.6*length(condition),pointsize = 1)
ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = length(names(res.ls)))
dev.off()
# Now check your DESeq2 folder. The plot should be in there now.

## Now we do the same for resLfs.ls results (shrunk Log2FC values)
gg <- list()
for(i in names(resLfs.ls)) {
  # keep the LFcut threshold from non shrunk Log2FC values
  LFcut = quantile(abs(res.ls[[i]][["log2FoldChange"]]), .9)
  res = resLfs.ls[[i]] # picking the i result table
  gg[[paste0("MA.",i)]] <- MAfun(res, Symbol = T, title = i, topgenes = 6, LFcut = LFcut, shrink = T)
  gg[[paste0("VO.",i)]] <- VulcFun(res, Symbol = T, title = i, topgenes = 6, LFcut = LFcut, shrink = T)
}
# png export of plot
png(paste0(substance,"_MA_Volcano_lfs.png"),units = "cm", bg = "white", res = 500,
    width = 7.33333*3, height = 7.6*length(condition),pointsize = 1)
ggpubr::ggarrange(plotlist = gg, ncol = 2, nrow = length(names(resLfs.ls)))
dev.off()

rm(gg,i,res,LFcut) # cleanup environment a little bit

# Q: Now compare the MA / Volcano plots between the "normal" and log2FC shrunk results.
#    How do the results differ?
############################



##############     NEW LINES OF CODE FROM HERE     ##################

### Venn plot ### --------------------------------------------------------------
# Based on what we have seen so far our data looks pretty good and reliable.
# Let's check  the numbers of differential expressed genes (DEGs) in the high 
# and low exposure condition again.

# list to store DEGs in
deg.ls <- lapply(res.ls, function(x){subset(x, padj <= p)})
# Here we are telling R to subset the df for all values which have a padj < p. p was defined as 0.05 previously
# Every gene in the tables stored in deg.ls has a padj <= 0.05

# Now let's check the number of DEGs in the highexposure condition
nrow(deg.ls[[grep("[Hh]igh",names(deg.ls))]])
# Wow! almost 5000 genes in here!!! 
# Q: Does this number make sense to be visualized in a heatmap? 
#    Or do you think it makes sense to visualize anything with ~ 5000 rows? 

# So in our case we might be only interested in genes that are found differentially
# regulated in both (HE & LE) conditions. 
# Q: Why would we be more interested in the common set of genes? 

# There is a really cool function in called ?intersect()
# Q: What does intersect() do? What can be union() and setdiff() be used for?
int <- intersect(row.names(deg.ls$LowExposure),
                 row.names(deg.ls$HighExposure))
# Q: Check the int object we just created. What does it contain?
# Q: what is the length(int)?

# Attention spoiler alert! Yes, the int contains all the gene IDs that were found
# dif.exp. in the common subset of high and low exposure condition. 
# We cann use a venn diagram to nicely visualize that:

# plot
myVenn <- function(deg.ls, shape = "circle", title = "") {
  col = topo.colors(length(names(deg.ls))) # color settings
  ven.ls <- lapply(deg.ls, rownames) # here we select only the rownames / Gene IDs from the tables
  venn <- eulerr::euler(ven.ls, shape = shape)
  s <- round(venn$stress,3)
  e <- round(venn$diagError,3)#
  plot(venn,
       fills = list(fill = col, alpha = .6), #labels = list(col = "black", font = 4),
       legend = list(col = "black", font = 4),
       main = paste0(title," - DEGs [Stress:",s," ; Diag.Er:",e,"]"),
       quantities = TRUE, shape = shape, lty = 1) #lty=0 for transparent
}
myVenn(deg.ls, title = substance)
# Great! Looking at the venn we see quite a big overlap between DEGs in HE and LE condition
# But let's be honest, ~450 genes is still too much for a nice heatmap.
# So let us filter our DEGs further. 

# The next thing we can filter for might be biological effect size (log2FC values)
# We can determine an arbitrary value, or we can compute one base on the LFcut distribution.
# To inspect the distr. of log2FC values you can plot them with:
par(mfrow=c(2,1))
hist(res.ls$LowExposure$log2FoldChange, breaks = 150, xlim = c(-5,5))
hist(res.ls$HighExposure$log2FoldChange, breaks = 150, xlim = c(-5,5))
# Q: What value are the log2FC values are centered around? 
# Q: What kind of data distribution do we have here?
# Q: How do the HE and LE log2FC values differ in their distribution?

# I have experienced in the past that using the upper 90% quantile (= top 10%) of genes with the absolute largest log2FC values)
# makes a great effect size cut off which scales with the "wideness" / variance of your log2FC distribution.
# We already know how to compute the quantile values. To get the top 90% we run:
LFcut.LE <- quantile(abs(res.ls$LowExposure$log2FoldChange), 0.9)
LFcut.LE
LFcut.HE <- quantile(abs(res.ls$HighExposure$log2FoldChange), 0.9)
# Q: why are we using the abs() expression here? 
# The LFcut.HE is actually quite large (> 1). To avoid the removal of precious datapoints
# we will manually set this value to 1 for now. But in fact there is NO textbook value. 
LFcut.HE = 1

# To filter out DEGs above or below the computed LFcut threshold run:
degFCcut.ls <- lapply(res.ls, function(X){
  lfcut = quantile(abs(X$log2FoldChange), 0.9)
  if(lfcut > 1){lfcut = 1} # in case lfcut greater 1, set to 1
  message(paste("Log2FC cut off:\t\t",round(lfcut,2)))
  x = subset(X, padj <= p)
  subset(x, abs(log2FoldChange) >= lfcut)
})
# Q: Inspect the list object degFCcut.ls. How many DEGs are left now in each condition?


# Let's plot the common set of DEGs after the LFcut filtering:
myVenn(degFCcut.ls, title = paste(coldata$Substance[1],"- LFcut"))
# This looks pretty good doesn't it? And in fact a number of 150 DEGs in the common set
# is much more reasonable to plot in the final heatmap.

# Extract these 150 gene IDs with:
tmp = lapply(degFCcut.ls, row.names)
names(tmp)
str(tmp) # Now we have a nice tidy list object with all ENSEMBL IDs of the DEGs for each condition
# To get the common overlap we simply run:
int <- Reduce(intersect, tmp)

# To do just one last checkup, let us have a look at the mean expression values of these 150 genes
tmp = apply(normMtx[int,],1,mean)
par(mfrow=c(1,1))
barplot(sort(tmp, decreasing = T), ylim = c(0,3000), main = "Mean expression of DEG selection")
abline(h=250, lwd=1.5, lty=2, col = "firebrick")
abline(h=500, lwd=1.5, lty=2, col = "blue")
# The two dashed lines indicate an arbitrary threshold we could set at i.e. 500 (blue)
# or 250 (red). So let us subset our final DEG selection one last time for genes that have
# a mean expression value of at least 100.
common <- names(tmp) # the 150 DEGs of the common set

minExp <- 200 # minimum baseMean gene count; You can play around with this parameter.
select <- names(tmp[which(tmp > minExp)]) # Final selection with mean expression cut off
select
# Yuhuuuuu! Here we have our final selection of 78 genes that:
# - were found dif.expressed in HE and LE condition
# - have a log2FC value greater than the respective cut off
# - a mean expression > 500

# Now finally! Let's make a heatmap ...
#################


### Plotting a heatmap ### -----------------------------------------------------------
# Here we will use the pheatmap function from the pheatmap package to plot the heatmap.
# There are many other tools out there. Feel free to explore :) 
?pheatmap

# very simple heatmap - only counts
pheatmap(countMtx[select,])

# now with annotation and gene symbol 
anno = subset(res.ls[[1]][select,], select=external_gene_name)
tmp = countMtx
tmp = merge(tmp,anno, by=0)
row.names(tmp) <- tmp$external_gene_name
tmp = subset(tmp, select = grep("R[1-9]",colnames(tmp)))
pheatmap(tmp, annotation_col = coldata[,c("Condition","Tank")])
# Well this looks sort of like a heatmap but probably not exactly what you expected.
# Q: What do the colors correspond here? 

# Plotting the vst transformed mean counts 
pheatmap(vst[select,], cluster_cols = T,
         annotation_col = coldata[,c("Condition","Tank")])
# Q: This looks nicer doesn't it? But can you tell from looking at the heatmap which genes are
# up and down regulated? What do the colors correspond to?
# You see, only having fancy colors doesn't make a plot informative. 

# Let us try to scale the rows.
# Basically what we are doing now is to compute the mean of each row then center (substract) every value
# of that row around it and then divide by the row's SD. This approach is called z-score transformation.
# Fancy name - simple math behind it ;).
pheatmap(vst[select,], scale = "row",
         annotation_col = coldata[,c("Condition","Tank")])
# Cool we are getting there , right?!
# Q: Can you tell now which genes are up and down regulated with respect to the control?

# Well, if you ask me this is not good enough yet. What about a heatmap where the values are
# centered not around the row's mean but around the control groups mean. 
# That way, the values would be ~ 0 and become negative when the expression is lower or positive
# when the expression is higher compared to the control. 
# Here is a function that will help you do that:
centCtrMean <- function(mtx, coldata){
  coldat <- coldata[colnames(mtx),]
  id <- rownames(coldat[grep("[Cc]ontrol", coldat$Condition),])
  ctrM <- apply(mtx[,id], 1, FUN = mean) #calcs the mean of control
  Sd <- apply(mtx, 1, FUN = sd) # calcs Sd of each row / gene
  (mtx - ctrM)/Sd # centers for the mean of control and scales for overall Sd of the row / gene
}

hMtx <- centCtrMean(vst[select,], coldata) #mean centered mtx with base mean cutoff for heatmap
hMtx.all <- centCtrMean(vst[common,], coldata) #mean centered mtx for heatmap without baseMean cutoff
# Now let's plot that:
pheatmap(hMtx)
# Q: Looking at the heatmap can you tell now which genes are up and down regulated? 

# Btw if you want no clustering of the columns you can run
pheatmap(hMtx, cluster_cols = F)
# to rearange the order of your columns try:
id.order <- c()
for(i in levels(as.factor(coldata$Substance))){
  for(k in levels(as.factor(coldata$Condition))){
    id <- rownames(coldata[coldata$Condition %in% k & coldata$Substance %in% i,])
    id.order <- append(id.order,id)
  }
}
hMtx <- hMtx[,id.order]
pheatmap(hMtx, cluster_cols = F, gaps_col = seq(1,ncol(hMtx)-3, by=3)+2,
         annotation_col = coldata[,c('Tank','Condition')])
# I think this is already a quite nice heatmap. All it needs is a little bit more fine tuning
# of the colors and some propper sample annotation. 

# Luckily I made a function for you ;):
myheat2 <- function(df, coldata, colclust = F, rowclust = T, Symbol = F,
                    title = "", distM = "euclidean", clustM = "average"){
  # Set anno colors
  col.Cond = colorRampPalette(RColorBrewer::brewer.pal(n=8, name="YlOrRd"))(length(levels(coldata$Condition)))
  col.Cond[1] <- "#56B4E9" #Defines color for control
  col.Tank = colorRampPalette(c("gray95","gray50"))(length(levels(coldata$Tank)))
  ann_colors = list(Tank = setNames(col.Tank, levels(coldata$Tank)),
                    Condition = setNames(col.Cond, levels(coldata$Condition)))
  # column gaps
  colGaps <- c()
  Conditions <- levels(coldata$Condition)
  for(k in Conditions[1:(length(Conditions)-1)]){
    x <- which(coldata[id.order,'Condition'] %in% k)
    x <- tail(x,1)
    colGaps <- append(colGaps,x)
  }
  # cluster fun
  callback <- function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  # set anno colors & breaks
  x <- 10 #length of colors above and below zero values
  x <- 2*x + 1
  col <- colorRampPalette(c("mediumblue","white","red2"))(x) #defines color palette in heatmap
  s <- sd(df[,rownames(coldata[which(coldata$Condition %in% condition[1]),])]) # sets sd around 0 values from control to white
  m <- mean(df[,rownames(coldata[which(coldata$Condition %in% condition[1]),])])
  myBreaks <- c(seq(min(df), m-s, length.out=ceiling(x/2)), 
                seq(m+s, max(df), length.out=floor(x/2)))
  # plot
  if(Symbol == T) {
    tmp = merge(df, res.ls[[1]], by = 0)[,-1]
    # if any external gene symbol is duplicated, order input by baseMean and
    # remove the lower expressed duplicate.
    if(any(duplicated(tmp$external_gene_name))) {
      tmp = dplyr::arrange(tmp, desc(baseMean))
      dup = which(duplicated(tmp$external_gene_name))
      dup = unique(tmp$external_gene_name[dup])
      tmp = tmp[!duplicated(tmp$external_gene_name), ]
      warning(paste("\nGene symbol",dup,"was duplicated. Lower expressed symbol was removed!"))
      row.names(tmp) <- tmp$external_gene_name
      } else { row.names(tmp) <- tmp$external_gene_name }
    df1 <- tmp[,1:ncol(df)]
    } else { df1 <- df }
  
  pheatmap::pheatmap(
    df1, angle_col = "90", gaps_col = if(colclust == T){NULL}else{colGaps},
    clustering_callback = callback, treeheight_col = 30, #default 50
    cluster_rows= rowclust, cluster_cols= colclust,
    clustering_distance_rows = distM, clustering_distance_cols = distM, 
    clustering_method = clustM,
    main = paste(title,distM, clustM,nrow(df1),"Genes"),
    show_rownames= T, show_colnames = T, breaks = myBreaks,
    annotation_col = coldata[,c('Tank','Condition')], annotation_colors = ann_colors, color = col
  )
}

# Here are your final heatmaps! 
#
myheat2(hMtx, coldata, Symbol = T, title = paste("min.av. gene count >",minExp,"-"))
myheat2(hMtx, coldata, Symbol = T, colclust = T, clustM = "average")

x <- hMtx.all[,id.order]
myheat2(x, coldata, Symbol = T, title = "Common set of DEGs")
myheat2(x, coldata, Symbol = T, colclust = T, clustM = "average")

# You can export them through:
pdf(file = "Heatmaps_commonDEGs.pdf", width = 5, height = 11, 
    onefile = T, bg = "transparent")
myheat2(hMtx, coldata, Symbol = T)
myheat2(hMtx, coldata, Symbol = T, colclust = T)
dev.off()
##########################

## Session Information & save Rdata  ## ------------
sink(paste0(home,"/DESeq2_SessionInfo.txt"))
print(date())
print(devtools::session_info())
sink()
save.image(paste0(home,"/DESeq2_pairwise.RData"))
###    END OF SCRIPT   ####