### DESeq2 in a nutshell ###
# Import
dds2 = DESeqDataSetFromMatrix(countData = countMtx, # <--- our count matrix
                              colData   = coldata,    # <--- our coldata file
                              design    = ~ Tank + Condition) #
# Filtering
dds2 = dds[rowSums(counts(dds2)) >= ncol(counts(dds2))+2, ]

# Run DESeq (normalisation, outlier detection, statistical testing & pvalue correction)
dds2 = DESeq(dds2, test = "Wald") 

# Inspect output
resultsNames(dds2)
# To get the 4th & 5th object simply type
resultsNames(dds2)[4:5]

# High exposure vs Control
res.HE = results(dds, name = resultsNames(dds2)[5], filterFun = ihw, alpha = p)
summary(res.HE)
head(res.HE)

# Low exposure vs Control
res.LE = results(dds, name = resultsNames(dds2)[4], filterFun = ihw, alpha = p)
summary(res.LE)
head(res.LE)