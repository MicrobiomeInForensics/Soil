## ----libraries--------------------------------------------------------------------------------------------------------------
library(decontam)
library(ggplot2)
library(phyloseq)


## ----phyloseq_object--------------------------------------------------------------------------------------------------------
#preparing the data
tab <- t(tab)
tax <- as.matrix(tax)
# checing that the data is in the correct configuration
identical(row.names(tab), row.names(tax))
identical(colnames(tab), row.names(metaData))

# Making the phyloseq object
ps <- phyloseq(otu_table(tab, taxa_are_rows=T), 
                  sample_data(metaData), 
                  tax_table(tax))
ps
head(sample_data(ps))


## ---------------------------------------------------------------------------------------------------------------------------
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type, shape = batch)) + geom_point()



## ---------------------------------------------------------------------------------------------------------------------------

contamdf.prev <- isContaminant(ps, method="prevalence", neg="ctrl", batch="batch")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))


contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="ctrl", threshold=0.5, batch="batch")
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant))

# Decided to use the stricter threshold (0.5-removes ASVs with 50% or more probability of being contaminant)


## ---------------------------------------------------------------------------------------------------------------------------
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$ctrl == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$ctrl == "FALSE", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") ## one sample is acting wierd? explanation of plot

## ---------------------------------------------------------------------------------------------------------------------------
ps_clean <- prune_taxa(!contamdf.prev05$contaminant, ps)
tab <- as.data.frame(t(as(otu_table(ps_clean), "matrix")))
dim(tab)
tax <- as.data.frame(as(tax_table(ps_clean), "matrix"))
dim(tax)
metaData <- as.data.frame(as(sample_data(ps_clean), "matrix"))

