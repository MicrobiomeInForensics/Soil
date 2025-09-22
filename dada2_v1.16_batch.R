#dada2 versjon kan kjores i R4.0 >
library('dada2')
packageVersion('dada2')
path<-'your_path_to_seqfiles'
list.files(path)
fnFs<-sort(list.files(path,pattern='_R1.fastq.gz',full.names=T))
fnRs<-sort(list.files(path,pattern='_R2.fastq.gz',full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_R1.fastq.gz"), `[`, 1)
pdf('./graph.pdf')
plotQualityProfile(fnFs[1:4])
#on cluster you can name the output plot: pdf('./QualProfF.pdf')
#on cluster remember to run dev.off() after commant for plot to be made
#truncate to 240bp
plotQualityProfile(fnRs[1:4])
#truncate to 160bp
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), trimLeft=c(19,20),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#~one min. for 4 samples PE on cluster with multithread=T
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
#<2 min. on cluster with multithread=T
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dev.off()

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#<30 sec. on cluster#only 1 unique seq for sample 3
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[3]]
dadaRs[[2]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)
write.table(track, "/cluster/projects/nn5017k/microbiome/Host_2024_prover/all_samples/track_sequence_loss_table.txt", sep = "\t")


taxa <- assignTaxonomy(seqtab.nochim, "/cluster/projects/nn5017k/microbiome/Host_2024_prover/A2_seqfiles/silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/cluster/projects/nn5017k/microbiome/Host_2024_prover/A2_seqfiles/silva/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
print(taxa.print)

#save.image(file='./dada2_v1.16.RData')

#write.table(seqtab.nochim,'./dada2_tab.txt')
write.table(taxa.print,'./dada2_taxonomy.txt',row.names=F)
seqs.out<-as.list(colnames(seqtab.nochim))
length(seqs.out)
names(seqs.out)<-paste('ASV',1:length(seqs.out),sep='')
library('seqinr')
write.fasta(seqs.out,names(seqs.out),'./ASV.fas',nbchar=500,as.string=T)

tab_out<-seqtab.nochim
colnames(tab_out)<-paste('ASV',1:length(seqs.out),sep='')
write.table(tab_out,'./dada2_tab.txt')

